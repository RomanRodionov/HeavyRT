#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>

#include "BVH2CommonRTStacklessLBVH.h"

void BVH2CommonRTStacklessLBVH::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
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

    if (v >= -1e-6f && u >= -1e-6f && (u + v <= 1.0f + 1e-6f) && t > tNear && t < pHit->t) // if (v > -1e-6f && u > -1e-6f && (u + v < 1.0f+1e-6f) && t > tMin && t < hit.t)
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

//extern bool g_debugPrint;

void BVH2CommonRTStacklessLBVH::BVH2TraverseC32(const float3 ray_pos, const float3 ray_dir, float tNear, 
                                                uint32_t instId, uint32_t geomId, uint32_t stack[STACK_SIZE], bool stopOnFirstHit,
                                                CRT_Hit *pHit)
{
  const uint32_t bvhOffset = m_bvhOffsets[geomId]; //

  int top = 0;
  uint32_t leftNodeOffset = 0; // + 1 because we skeep root node and directly go to it's children
  #ifdef ENABLE_METRICS
  uint32_t leftNodeOffsetOld = leftNodeOffset;
  #endif

  const float3 rayDirInv = SafeInverse(ray_dir);

  while (top >= 0 && !(stopOnFirstHit && pHit->primId != uint32_t(-1)))
  {
    //if(g_debugPrint)
    //  std::cout << "leftNodeOffset = " << leftNodeOffset << std::endl;

    bool needStackPop = (leftNodeOffset == 0xFFFFFFFF); // for empty leaves
    if ((leftNodeOffset & LEAF_BIT) == 0)
    {
      #ifdef ENABLE_METRICS
      m_stats.NC  += 2;
      m_stats.BLB += 2*sizeof(BVHNode);
      if (leftNodeOffsetOld != leftNodeOffset) {
        for (int i = 0; i < TREELET_ARR_SIZE; i++) {
          if (std::abs(double(leftNodeOffset) - double(leftNodeOffsetOld)) * double(sizeof(BVHNode)) >= (double)treelet_sizes[i])
            m_stats.LJC[i]++;
          const uint32_t oldCacheLineId = uint32_t(leftNodeOffsetOld * sizeof(BVHNode)) / uint32_t(treelet_sizes[i]);
          const uint32_t newCacheLineId = uint32_t(leftNodeOffset * sizeof(BVHNode)) / uint32_t(treelet_sizes[i]);
          if (oldCacheLineId != newCacheLineId){
            m_stats.CMC[i]++;
            m_stats.WSS[i].insert(newCacheLineId);
          }
        }
        leftNodeOffsetOld = leftNodeOffset;
      }
      #endif
      
      //auto nodeAddrB = reinterpret_cast<uintptr_t>(m_allNodes.data() + bvhOffset);
      //auto nodeAddr  = reinterpret_cast<uintptr_t>(m_allNodes.data() + bvhOffset + leftNodeOffset);
      //assert(nodeAddrB % 64 == 0);
      //assert(nodeAddr  % 64 == 0);

      const BVHNode node0 = m_allNodes[bvhOffset + leftNodeOffset + 0];
      const BVHNode node1 = m_allNodes[bvhOffset + leftNodeOffset + 1];

      const float2 tm0 = RayBoxIntersection2(ray_pos, rayDirInv, node0.boxMin, node0.boxMax);
      const float2 tm1 = RayBoxIntersection2(ray_pos, rayDirInv, node1.boxMin, node1.boxMax);

      const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= tNear) && (tm0.x <= pHit->t);
      const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= tNear) && (tm1.x <= pHit->t);
      needStackPop         = (!hitChild0 && !hitChild1);

      //if(g_debugPrint)
      //{
      //  std::cout << "node0.boxMin = (" << node0.boxMin.x << ", " << node0.boxMin.y << ", " << node0.boxMin.z << ")" << std::endl;
      //  std::cout << "node0.boxMax = (" << node0.boxMax.x << ", " << node0.boxMax.y << ", " << node0.boxMax.z << ")" << std::endl;
      //}

      // traversal decision
      //
      leftNodeOffset = hitChild0 ? node0.leftOffset : node1.leftOffset;

      if (hitChild0 && hitChild1)
      {
        leftNodeOffset = (tm0.x <= tm1.x) ? node0.leftOffset : node1.leftOffset; // GPU style branch
        stack[top]     = (tm0.x <= tm1.x) ? node1.leftOffset : node0.leftOffset; // GPU style branch
        top++;
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
      }
    } // end while (searchingForLeaf)

    // leaf node, intersect triangles
    //
    if ((leftNodeOffset & LEAF_BIT) != 0 && leftNodeOffset != 0xFFFFFFFF)
    {
      //if(g_debugPrint)
      //  std::cout << "seek for intersection at " << leftNodeOffset << std::endl;
      const uint32_t start = EXTRACT_START(leftNodeOffset);
      const uint32_t count = EXTRACT_COUNT(leftNodeOffset);
      IntersectAllPrimitivesInLeaf(ray_pos, ray_dir, tNear, instId, geomId, start, count, pHit);
      needStackPop = true;
      #ifdef ENABLE_METRICS
      m_stats.LC++;
      m_stats.LC2++;
      m_stats.TC+=count;
      m_stats.BLB+=count*(3*sizeof(uint32_t) + 3*sizeof(float4));
      #endif
    }

    // continue BVH traversal
    //
    if(needStackPop)
    {
      top--;
      leftNodeOffset = stack[std::max(top,0)];
    }

    #ifdef ENABLE_METRICS
    m_stats.SOC++;
    m_stats.SBL+=sizeof(uint32_t); 
    #endif
  } // end while (top >= 0)
}
             
CRT_Hit BVH2CommonRTStacklessLBVH::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  const bool stopOnFirstHit = (dirAndFar.w <= 0.0f);
  if(stopOnFirstHit)
    dirAndFar.w *= -1.0f;

  #ifdef ENABLE_METRICS
  ResetVarLC();
  #endif

  uint32_t stack[STACK_SIZE];

  // (2) process all intersected BLAS
  //
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
  
      BVH2TraverseC32(ray_pos, ray_dir, posAndNear.w, instId, geomId, stack, stopOnFirstHit,
                      &hit);
    }

  } while (nodeIdx < 0xFFFFFFFE && !(stopOnFirstHit && hit.primId != uint32_t(-1)));

  #ifdef ENABLE_METRICS
  m_stats.raysNumber++;
  #endif
  
  #ifdef REMAP_PRIM_ID
  if(hit.geomId < uint32_t(-1)) 
  {
    const uint2 geomOffsets = m_geomOffsets[hit.geomId];
    hit.primId = m_primIndices[geomOffsets.x/3 + hit.primId];
  }
  #endif

  return hit;
}

bool BVH2CommonRTStacklessLBVH::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  dirAndFar.w *= -1.0f;
  CRT_Hit hit = RayQuery_NearestHit(posAndNear, dirAndFar);
  return (hit.geomId != uint32_t(-1));
}
