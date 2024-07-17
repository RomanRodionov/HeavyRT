#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <chrono>

#include "BVH2FatRT.h"

void BVH2FatRT::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
                                             float tNear, uint32_t instId, uint32_t geomId,
                                             uint32_t a_start, uint32_t a_count,
                                             CRT_Hit *pHit)
{
  const uint2 a_geomOffsets = m_geomOffsets[geomId];

  for (uint32_t triId = a_start; triId < a_start + a_count; triId++)
  {
    const uint32_t A = m_indices[a_geomOffsets.x + triId*3 + 0];
    const uint32_t B = m_indices[a_geomOffsets.x + triId*3 + 1];
    const uint32_t C = m_indices[a_geomOffsets.x + triId*3 + 2];

    const float3 A_pos = to_float3(m_vertPos[a_geomOffsets.y + A]);
    const float3 B_pos = to_float3(m_vertPos[a_geomOffsets.y + B]);
    const float3 C_pos = to_float3(m_vertPos[a_geomOffsets.y + C]);

    const float3 edge1 = B_pos - A_pos;
    const float3 edge2 = C_pos - A_pos;
    const float3 pvec = cross(ray_dir, edge2);
    const float3 tvec = ray_pos - A_pos;
    const float3 qvec = cross(tvec, edge1);

    const float invDet = 1.0f / dot(edge1, pvec);
    const float v = dot(tvec, pvec) * invDet;
    const float u = dot(qvec, ray_dir) * invDet;
    const float t = dot(edge2, qvec) * invDet;

    if (v >= -1e-6f && u >= -1e-6f && (u + v <= 1.0f + 1e-6f) && t > tNear && t < pHit->t) 
    {
      pHit->t = t;
      pHit->primId = triId;
      pHit->instId = instId;
      pHit->geomId = geomId;
      pHit->coords[0] = u;
      pHit->coords[1] = v;
    }
  }
}

void BVH2FatRT::BVH2TraverseF32(const float3 ray_pos, const float3 ray_dir, float tNear,
                                uint32_t instId, uint32_t geomId, uint32_t stack[STACK_SIZE], bool stopOnFirstHit,
                                CRT_Hit* pHit)
{
  const uint32_t bvhOffset = m_bvhOffsets[geomId];

  int top = 0;
  uint32_t leftNodeOffset = 0;

#ifdef ENABLE_METRICS
  uint32_t leftNodeOffsetOld = leftNodeOffset;
#endif

  const float3 rayDirInv = SafeInverse(ray_dir);
  while (top >= 0 && !(stopOnFirstHit && pHit->primId != uint32_t(-1)))
  {
    while (top >= 0 && ((leftNodeOffset & LEAF_BIT) == 0))
    {
      #ifdef ENABLE_METRICS
      m_stats.NC  += 2;
      m_stats.BLB += 1 * sizeof(BVHNodeFat);
      if (leftNodeOffsetOld != leftNodeOffset) {
        for (int i = 0; i < TREELET_ARR_SIZE; i++) {
          if (std::abs(double(leftNodeOffset) - double(leftNodeOffsetOld)) * double(sizeof(BVHNodeFat)) >= (double)treelet_sizes[i])
            m_stats.LJC[i]++;
          const uint32_t oldCacheLineId = uint32_t(leftNodeOffsetOld * sizeof(BVHNodeFat)) / uint32_t(treelet_sizes[i]);
          const uint32_t newCacheLineId = uint32_t(leftNodeOffset * sizeof(BVHNodeFat)) / uint32_t(treelet_sizes[i]);
          if (oldCacheLineId != newCacheLineId){
            m_stats.CMC[i]++;
            m_stats.WSS[i].insert(newCacheLineId);
          }
        }
        leftNodeOffsetOld = leftNodeOffset;
      }
      #endif

      const BVHNodeFat fatNode = m_allNodesFat[bvhOffset + leftNodeOffset];

      const float3 leftBoxMin  = to_float3(fatNode.lmin_xyz_rmax_x);
      const float3 leftBoxMax  = to_float3(fatNode.lmax_xyz_rmax_y);
      const float3 rightBoxMin = to_float3(fatNode.rmin_xyz_rmax_z);
      const float3 rightBoxMax = float3(fatNode.lmin_xyz_rmax_x.w, fatNode.lmax_xyz_rmax_y.w, fatNode.rmin_xyz_rmax_z.w);

      const uint32_t node0_leftOffset = fatNode.offs_left;
      const uint32_t node1_leftOffset = fatNode.offs_right;

      const float2 tm0 = RayBoxIntersection2(ray_pos, rayDirInv, leftBoxMin, leftBoxMax);
      const float2 tm1 = RayBoxIntersection2(ray_pos, rayDirInv, rightBoxMin, rightBoxMax);

      const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= tNear) && (tm0.x <= pHit->t);
      const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= tNear) && (tm1.x <= pHit->t);

      // traversal decision
      //
      leftNodeOffset = hitChild0 ? node0_leftOffset : node1_leftOffset;

      if (hitChild0 && hitChild1)
      {
        leftNodeOffset = (tm0.x <= tm1.x) ? node0_leftOffset : node1_leftOffset; // GPU style branch
        stack[top]     = (tm0.x <= tm1.x) ? node1_leftOffset : node0_leftOffset; // GPU style branch
        top++;
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL += sizeof(uint32_t);
        #endif
      }

      if (!hitChild0 && !hitChild1) // both miss, stack.pop()
      {
        top--;
        leftNodeOffset = stack[std::max(top,0)];
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL += sizeof(uint32_t);
        #endif
      }

    } // end while (searchingForLeaf)

    // leaf node, intersect triangles
    //
    if (top >= 0 && leftNodeOffset != 0xFFFFFFFF)
    {
      const uint32_t start = EXTRACT_START(leftNodeOffset);
      const uint32_t count = EXTRACT_COUNT(leftNodeOffset);
      IntersectAllPrimitivesInLeaf(ray_pos, ray_dir, tNear, instId, geomId, start, count, pHit);
      #ifdef ENABLE_METRICS
      m_stats.LC++;
      m_stats.LC2++;
      m_stats.TC += count;
      m_stats.BLB += count * (3 * sizeof(uint32_t) + 3 * sizeof(float4));
      #endif
    }

    // continue BVH traversal
    //
    top--;
    leftNodeOffset = stack[std::max(top,0)];

    #ifdef ENABLE_METRICS
    m_stats.SOC++;
    m_stats.SBL += sizeof(uint32_t);
    #endif
  } // end while (top >= 0)

}

CRT_Hit BVH2FatRT::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  bool stopOnFirstHit = (dirAndFar.w <= 0.0f);
  if(stopOnFirstHit)
    dirAndFar.w *= -1.0f;

  #ifdef ENABLE_METRICS
  ResetVarLC();
  m_stats.raysNumber++;
  #endif

  uint32_t stack[STACK_SIZE];

  CRT_Hit hit;
  hit.t      = dirAndFar.w;
  hit.primId = uint32_t(-1);
  hit.instId = uint32_t(-1);
  hit.geomId = uint32_t(-1);

  const float3 rayDirInv = SafeInverse(to_float3(dirAndFar));
  uint32_t nodeIdx = 0;
  do
  {
    uint32_t travFlags  = 0;
    uint32_t leftOffset = 0;
    do
    {
      const BVHNode currNode = m_nodesTLAS[nodeIdx];
      const float2 boxHit    = RayBoxIntersection2(to_float3(posAndNear), rayDirInv, currNode.boxMin, currNode.boxMax);
      const bool intersects  = (boxHit.x <= boxHit.y) && (boxHit.y > posAndNear.w) && (boxHit.x < hit.t); // (tmin <= tmax) && (tmax > 0.f) && (tmin < curr_t)
      
      #ifdef ENABLE_METRICS
      m_stats.NC  += 1;
      m_stats.BLB += 1 * sizeof(BVHNode);
      #endif

      travFlags  = (currNode.leftOffset & LEAF_BIT) | uint32_t(intersects); // travFlags  = (((currNode.leftOffset & LEAF_BIT) == 0) ? 0 : LEAF_BIT) | (intersects ? 1 : 0);
      leftOffset = currNode.leftOffset;
      nodeIdx    = isLeafOrNotIntersect(travFlags) ? currNode.escapeIndex : leftOffset;

    } while (notLeafAndIntersect(travFlags) && nodeIdx != 0 && nodeIdx < 0xFFFFFFFE); 
     
    if(isLeafAndIntersect(travFlags)) 
    {
      const uint32_t instId = EXTRACT_START(leftOffset);
      const uint32_t geomId = m_geomIdByInstId[instId];
  
      // transform ray with matrix to local space
      //
      const float3 ray_pos = matmul4x3(m_instMatricesInv[instId], to_float3(posAndNear));
      const float3 ray_dir = matmul3x3(m_instMatricesInv[instId], to_float3(dirAndFar)); // DON'T NORMALIZE IT !!!! When we transform to local space of node, ray_dir must be unnormalized!!!
  
      BVH2TraverseF32(ray_pos, ray_dir, posAndNear.w, instId, geomId, stack, stopOnFirstHit, &hit);
    }

  } while (nodeIdx < 0xFFFFFFFE && !(stopOnFirstHit && hit.primId != uint32_t(-1))); //

  #ifdef REMAP_PRIM_ID
  if(hit.geomId < uint32_t(-1)) 
  {
    const uint2 geomOffsets = m_geomOffsets[hit.geomId];
    hit.primId = m_primIndices[geomOffsets.x/3 + hit.primId];
  }
  #endif 
  
  return hit;
}

bool BVH2FatRT::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  dirAndFar.w *= -1.0f;
  CRT_Hit hit = RayQuery_NearestHit(posAndNear, dirAndFar);
  return (hit.geomId != uint32_t(-1));
}
