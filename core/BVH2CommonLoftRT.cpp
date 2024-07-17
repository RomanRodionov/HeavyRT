#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>

#include "BVH2CommonLoftRT.h"

void BVH2CommonLoftRT::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
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

CRT_Hit BVH2CommonLoftRT::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
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
  
  ////////////////////////////////////////////////////////////////////////////// instancing variables are a bit more than common
  int instTop  = 0;
  uint32_t instId = uint32_t(-1);
  
  float3 ray_pos = to_float3(posAndNear);
  float3 ray_dir = to_float3(dirAndFar);
  float3 inv_dir = SafeInverse(ray_dir);

  uint32_t bvhOffset = m_tlasOffset;
  ////////////////////////////////////////////////////////////////////////////// instancing variables are a bit more than common

  int top = 0;
  uint32_t leftNodeOffset = 0; 

  while (top >= 0 && !(stopOnFirstHit && hit.primId != uint32_t(-1)))
  {
    while (top >= 0 && ((leftNodeOffset & LEAF_BIT) == 0))
    {
      //if(g_debugPrint)
      //  std::cout << "leftNodeOffset = " << leftNodeOffset << std::endl;

      #ifdef ENABLE_METRICS
      m_stats.NC  += 2;
      m_stats.BLB += 2*sizeof(BVHNode);
      #endif
      
      //auto nodeAddrB = reinterpret_cast<uintptr_t>(m_allNodes.data() + bvhOffset);
      //auto nodeAddr  = reinterpret_cast<uintptr_t>(m_allNodes.data() + bvhOffset + leftNodeOffset);
      //assert(nodeAddrB % 64 == 0);
      //assert(nodeAddr  % 64 == 0);

      const BVHNode node0 = m_allNodes[bvhOffset + leftNodeOffset + 0];
      const BVHNode node1 = m_allNodes[bvhOffset + leftNodeOffset + 1];

      const float2 tm0 = RayBoxIntersection2(ray_pos, inv_dir, node0.boxMin, node0.boxMax);
      const float2 tm1 = RayBoxIntersection2(ray_pos, inv_dir, node1.boxMin, node1.boxMax);

      const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= posAndNear.w) && (tm0.x <= hit.t);
      const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= posAndNear.w) && (tm1.x <= hit.t);

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

      if (!hitChild0 && !hitChild1) // both miss, stack.pop()
      {
        top--;
        leftNodeOffset = stack[std::max(top,0)];
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
      }

      if (top >= 0 && top < instTop && bvhOffset != m_tlasOffset)
      {
        ray_pos   = to_float3(posAndNear);
        ray_dir   = to_float3(dirAndFar);
        inv_dir   = SafeInverse(ray_dir);
        bvhOffset = m_tlasOffset;
      }

    } // end while (searchingForLeaf)

    if (top >= 0 && leftNodeOffset != 0xFFFFFFFF && bvhOffset != m_tlasOffset)  // leaf node of BLAS, intersect triangles
    {
      const uint32_t start  = EXTRACT_START(leftNodeOffset);
      const uint32_t count  = EXTRACT_COUNT(leftNodeOffset);
      const uint32_t geomId = m_geomIdByInstId[instId];

      IntersectAllPrimitivesInLeaf(ray_pos, ray_dir, posAndNear.w, instId, geomId, start, count, &hit); 
      
      //if(g_debugPrint)
      //{
      //  std::cout << "seek for intersection at " << leftNodeOffset << std::endl;
      //  std::cout << "hit.t      = " << hit.t << std::endl;
      //  std::cout << "hit.primId = " << hit.primId << std::endl;
      //}

      #ifdef ENABLE_METRICS
      m_stats.LC++;
      m_stats.LC2++;
      m_stats.TC+=count;
      m_stats.BLB+=count*(3*sizeof(uint32_t) + 3*sizeof(float4));
      #endif
    }
    else if (top >= 0 && bvhOffset == m_tlasOffset)                            // leaf node of BLAS, intersect BLAS next
    {
      instId           = EXTRACT_START(leftNodeOffset);
      bvhOffset        = m_bvhOffsets[m_geomIdByInstId[instId]];
      
      //if(g_debugPrint)
      //  std::cout << "TLAS -> BLAS for inst(" << instId << ") at " << leftNodeOffset << std::endl;
        
      leftNodeOffset   = 0;

      ray_pos = matmul4x3(m_instMatricesInv[instId], to_float3(posAndNear));
      ray_dir = matmul3x3(m_instMatricesInv[instId], to_float3(dirAndFar));
      inv_dir = SafeInverse(ray_dir);
      instTop = top;
    }
    
    if(top >= 0 && (leftNodeOffset != 0 || leftNodeOffset == 0xFFFFFFFF)) // continue BVH traversal, except the cases when we have just enter from TLAS to BLAS (see upper code)
    {                                                                     // (leftNodeOffset == 0xFFFFFFFF); // for empty leaves
      top--;
      leftNodeOffset = stack[std::max(top,0)];

      if (top < instTop && bvhOffset != m_tlasOffset)
      {
        ray_pos   = to_float3(posAndNear);
        ray_dir   = to_float3(dirAndFar);
        inv_dir   = SafeInverse(ray_dir);
        bvhOffset = m_tlasOffset;
      }
      
      #ifdef ENABLE_METRICS
      m_stats.SOC++;
      m_stats.SBL+=sizeof(uint32_t); 
      #endif
    }
    
  } // end while (top >= 0)
  
  #ifdef REMAP_PRIM_ID
  if(hit.geomId < uint32_t(-1)) 
  {
    const uint2 geomOffsets = m_geomOffsets[hit.geomId];
    hit.primId = m_primIndices[geomOffsets.x/3 + hit.primId];
  }
  #endif

  return hit;
}

bool BVH2CommonLoftRT::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  dirAndFar.w *= -1.0f;
  CRT_Hit hit = RayQuery_NearestHit(posAndNear, dirAndFar);
  return (hit.geomId != uint32_t(-1));
}
