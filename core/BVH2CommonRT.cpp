#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>

#include "BVH2CommonRT.h"

void BVH2CommonRT::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
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

void BVH2CommonRT::BVH2TraverseC32(const float3 ray_pos, const float3 ray_dir, float tNear, 
                                   uint32_t instId, uint32_t geomId, uint32_t stack[STACK_SIZE], 
                                   CRT_Hit *pHit)
{
  const uint32_t bvhOffset = m_bvhOffsets[geomId];

  int top = 0;
  uint32_t leftNodeOffset = 0; // + 1 because we skeep root node and directly go to it's children
  #ifdef ENABLE_METRICS
  uint32_t leftNodeOffsetOld = leftNodeOffset;
  #endif

  const float3 rayDirInv = SafeInverse(ray_dir);

  while (top >= 0)
  {
    bool needStackPop = false;
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
    } 
    else if ((leftNodeOffset & LEAF_BIT) != 0 && leftNodeOffset != 0xFFFFFFFF)  // leaf node, intersect triangles
    {
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
      #ifdef ENABLE_METRICS
      m_stats.SOC++;
      m_stats.SBL+=sizeof(uint32_t); 
      #endif
    }
  } // end while (top >= 0)
}
             
uint32_t BVH2CommonRT::LBVH2Traverse(float4 posAndNear, float4 dirAndFar, uint32_t stack[STACK_SIZE],
                                     BoxHit out_hits[LBVH_MAXHITS])
{
  if((m_nodesTLAS[0].leftOffset & LEAF_BIT) != 0) // if root is leaf
  {
    out_hits[0] = make_BoxHit(EXTRACT_START(m_nodesTLAS[0].leftOffset), 1e37f);
    return 1;
  }  
  
  const float3 rayDirInv = SafeInverse(to_float3(dirAndFar));
  const float  tNear     = posAndNear.w;
  const float  tFar      = dirAndFar.w;

  int   topRes  = 0;
  int   top     = 0;
  
  uint32_t leftNodeOffset  = m_nodesTLAS[0].leftOffset;  // assume we intersect root and the root is never leaf
  uint32_t rightNodeOffset = m_nodesTLAS[0].escapeIndex;
  
  // while-while traversal

  while(top >= 0) 
  {
    while(top >= 0 && ((leftNodeOffset & LEAF_BIT) == 0)) // searching for leaf
    {
      const BVHNode node0 = m_nodesTLAS[leftNodeOffset];
      const BVHNode node1 = m_nodesTLAS[rightNodeOffset];
      
      #ifdef ENABLE_METRICS
      m_stats.NC  += 2;
      m_stats.BLB += 2 * sizeof(BVHNode);
      #endif

      const float2 tm0 = RayBoxIntersection2(to_float3(posAndNear), rayDirInv, node0.boxMin, node0.boxMax);
      const float2 tm1 = RayBoxIntersection2(to_float3(posAndNear), rayDirInv, node1.boxMin, node1.boxMax);
      
      const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= tNear) && (tm0.x <= tFar);
      const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= tNear) && (tm1.x <= tFar);
      const bool leftFirst = (tm0.x <= tm1.x);

      // traversal decision
      //
      leftNodeOffset  = hitChild0 ? node0.leftOffset  : node1.leftOffset;
      rightNodeOffset = hitChild0 ? node0.escapeIndex : node1.escapeIndex;

      if (hitChild0 && hitChild1)
      {
        leftNodeOffset  = leftFirst ? node0.leftOffset  : node1.leftOffset;
        rightNodeOffset = leftFirst ? node0.escapeIndex : node1.escapeIndex;
        stack[top+2+0]  = leftFirst ? node1.leftOffset  : node0.leftOffset;
        stack[top+2+1]  = leftFirst ? node1.escapeIndex : node0.escapeIndex;
        top+=2;
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
      }

      if (!hitChild0 && !hitChild1) // both miss, stack.pop()
      {
        top-=2;
        leftNodeOffset  = stack[top+2+0];
        rightNodeOffset = stack[top+2+1];
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
      }
    } 

    while((leftNodeOffset & LEAF_BIT) != 0 && leftNodeOffset != 0xFFFFFFFF && top >= 0)
    {
      out_hits[topRes]= make_BoxHit(EXTRACT_START(leftNodeOffset), 0.0f); 
      topRes++;
      
      top-=2;
      leftNodeOffset  = stack[top+2+0];
      rightNodeOffset = stack[top+2+1];
    }

  } // travse TLAS again
  
  return topRes;
} 

CRT_Hit BVH2CommonRT::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  #ifdef ENABLE_METRICS
  ResetVarLC();
  #endif

  BoxHit boxMinHits[LBVH_MAXHITS];
  uint32_t stack[STACK_SIZE];

  // (1) process TLAS
  //
  uint32_t boxesNum = LBVH2Traverse(posAndNear, dirAndFar, stack,
                                    boxMinHits);

  // (2) process all intersected BLAS
  //
  CRT_Hit hit;
  hit.t      = dirAndFar.w;
  hit.primId = uint32_t(-1);
  hit.instId = uint32_t(-1);
  hit.geomId = uint32_t(-1);

  for (int boxId = int(boxesNum)-1; boxId >= 0; boxId--)
  {
    if (boxMinHits[boxId].tHit <= hit.t)
    {
      const uint32_t instId = boxMinHits[boxId].id;
      const uint32_t geomId = m_geomIdByInstId[instId];
  
      // transform ray with matrix to local space
      //
      const float3 ray_pos = matmul4x3(m_instMatricesInv[instId], to_float3(posAndNear));
      const float3 ray_dir = matmul3x3(m_instMatricesInv[instId], to_float3(dirAndFar)); // DON'T NORMALIZE IT !!!! When we transform to local space of node, ray_dir must be unnormalized!!!
  
      BVH2TraverseC32(ray_pos, ray_dir, posAndNear.w, instId, geomId, stack, &hit);
    }
  }

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

bool BVH2CommonRT::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  #ifdef ENABLE_METRICS
  m_stats.raysNumber++;
  #endif

  // (2) If any hit is found, immediately return true.
  std::cout << "[BVH2CommonRT::RayQuery_AnyHit]: "
            << "not implemeted!" << std::endl;

  return false;
}
