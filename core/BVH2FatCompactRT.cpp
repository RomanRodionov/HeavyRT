#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <chrono>

#include "BVH2FatCompactRT.h"

void BVH2FatRTCompact::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
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

static inline uint32_t GetTreeletRootOffset(uint32_t a_dataOffset, uint32_t a_log2_group_size)
{
  return 1 + ((a_dataOffset >> a_log2_group_size) << a_log2_group_size);
}

static constexpr uint MAX_OFFSET_MASK     = 0x07FFFFFF;  // max 27 bits for offset
static constexpr uint MAX_TRG_START_MASK  = 0x007FFFFF;  // max 23 bits for the first triangle index
static constexpr uint MAX_TRG_COUNT_MASK  = 0x0000000F;  // max 4 bits for the number of triangles in a leaf

static inline float2 unpackFloat2x16(uint32_t v)
{
  LiteMath::half xy[2];
  memcpy(xy, &v, sizeof(uint32_t));
  return float2(xy[0], xy[1]);
}

static constexpr uint32_t LEFT_FIRST = 1;
static constexpr uint32_t HIT_LEFT   = 2;
static constexpr uint32_t HIT_RIGHT  = 4;
static constexpr uint32_t HIT_BOTH   = HIT_LEFT | HIT_RIGHT;

static inline bool hitBoth  (uint32_t a_travFlags) { return ((a_travFlags & HIT_BOTH) == HIT_BOTH); }
static inline bool missBoth (uint32_t a_travFlags) { return ((a_travFlags & HIT_BOTH) == 0); }


static inline uint3 TraversalMath(float3 ray_pos, float3 rayDirInv, uint4 data, float3 bb_origin, float3 bb_ssize, float tNear, float tFar)
{
  float2 tm0;
  {
    float3 leftBoxMin;
    leftBoxMin.x  = bb_origin.x + float(  data.x & 0x0000003F      )  * bb_ssize.x;
    leftBoxMin.y  = bb_origin.y + float( (data.x & 0x00000FC0) >> 6)  * bb_ssize.y;
    leftBoxMin.z  = bb_origin.z + float( (data.x & 0x0003F000) >> 12) * bb_ssize.z; 
    
    float3 leftBoxMax;
    leftBoxMax.x  = bb_origin.x + float( (data.x & 0x00FC0000) >> 18) * bb_ssize.x;
    leftBoxMax.y  = bb_origin.y + float( (data.x & 0x3F000000) >> 24) * bb_ssize.y;
    leftBoxMax.z  = bb_origin.z + float(((data.x & 0xC0000000) >> 30) | ((data.y & 0x0000000F) << 2) ) * bb_ssize.z;
    
    tm0 = RayBoxIntersection2(ray_pos, rayDirInv, leftBoxMin, leftBoxMax);
  }

  float2 tm1;
  {
    float3 rightBoxMin;
    rightBoxMin.x = bb_origin.x + float( (data.y & 0x000003F0) >> 4 ) * bb_ssize.x; 
    rightBoxMin.y = bb_origin.y + float( (data.y & 0x0000FC00) >> 10) * bb_ssize.y;
    rightBoxMin.z = bb_origin.z + float( (data.y & 0x003F0000) >> 16) * bb_ssize.z;
    
    float3 rightBoxMax;
    rightBoxMax.x = bb_origin.x + float( (data.y & 0x0FC00000) >> 22) * bb_ssize.x;
    rightBoxMax.y = bb_origin.y + float(((data.y & 0xF0000000) >> 28) | ((data.z & 0x00000003) << 4)) * bb_ssize.y;
    rightBoxMax.z = bb_origin.z + float( (data.z & 0x000000FC) >> 2 ) * bb_ssize.z;

    tm1 = RayBoxIntersection2(ray_pos, rayDirInv, rightBoxMin, rightBoxMax);
  }

  const bool hitChild0     = (tm0.x <= tm0.y) && (tm0.y >= tNear) && (tm0.x <= tFar);
  const bool hitChild1     = (tm1.x <= tm1.y) && (tm1.y >= tNear) && (tm1.x <= tFar);
  const bool hit_first     = (tm0.x <= tm1.x);
  const uint32_t travFlags = uint(hit_first) | (uint(hitChild0) << 1) | (uint(hitChild1) << 2);

  uint32_t offs_left = 0;   // left part
  if ((data.w & 0x00000008) != 0)  // check leaf bit (data[3]: bit:3)
  {
    const uint start = (data.z >> 8) & MAX_TRG_START_MASK;
    const uint count = ((data.z >> 31) | (data.w << 1)) & MAX_TRG_COUNT_MASK;
    offs_left        = PackOffsetAndSize(start, count);
  }
  else
  {
    offs_left = ((data.z >> 8) | (data.w << 24)) & MAX_OFFSET_MASK;
  }
  
  uint32_t offs_right = 0;  // right part
  if ((data.w & 0x80000000) != 0)  // check leaf bit (data[3]: bit:31)
  {
    const uint start = (data.w >> 4) & MAX_TRG_START_MASK;
    const uint count = (data.w >> 27) & MAX_TRG_COUNT_MASK;
    offs_right       = PackOffsetAndSize(start, count);
  }
  else
    offs_right = (data.w >> 4) & MAX_OFFSET_MASK;

  return uint3(offs_left, offs_right, travFlags);
}


void BVH2FatRTCompact::BVH2TraverseF32(const float3 ray_pos, const float3 ray_dir, float tNear, 
                                       uint32_t instId, uint32_t geomId, uint32_t stack[STACK_SIZE], bool stopOnFirstHit,
                                       CRT_Hit *pHit)
{
  const uint32_t bvhOffset = m_bvhOffsets[geomId]; // m_allNodesFat16B.Offset(geomId);

  int top = 0;
  uint32_t leftNodeOffset = 1; // The first block(s) in array are occupied by the root's BB, so the first node is at '1'
  const float3 rayDirInv  = SafeInverse(ray_dir);

  uint32_t currTreeletRootOffset = 0;  // The current treelet root
  float3 trBoxMin, trBoxSize;    // The data on the treelet's bounding box

  while (top >= 0 && !(stopOnFirstHit && pHit->primId != uint32_t(-1)))
  {
    while (top >= 0 && ((leftNodeOffset & LEAF_BIT) == 0))
    {
      const uint32_t tr_rootOffset = GetTreeletRootOffset(leftNodeOffset, m_log2_group_size);

      #ifdef ENABLE_METRICS
      m_stats.NC += 2;
      if (currTreeletRootOffset != tr_rootOffset)
        m_stats.BLB += 2 * sizeof(uint4);
      else
        m_stats.BLB += 1 * sizeof(uint4); 
      #endif
      
      if (currTreeletRootOffset != tr_rootOffset)
      {
        currTreeletRootOffset = tr_rootOffset;
        const uint4 trRoot    = m_allNodesFat16B[bvhOffset + tr_rootOffset - 1];

        const float2 t0 = float2(unpackFloat2x16(trRoot.x));
        const float2 t1 = float2(unpackFloat2x16(trRoot.y));
        const float2 t2 = float2(unpackFloat2x16(trRoot.z));
        const float2 t3 = float2(unpackFloat2x16(trRoot.w));  
      
        trBoxMin.x = t0.x;
        trBoxMin.y = t0.y;
        trBoxMin.z = t1.x;
              
        trBoxSize.x = t2.x;
        trBoxSize.y = t2.y;
        trBoxSize.z = t3.x;
      }
      
      const uint3 leftOffsets = TraversalMath(ray_pos, rayDirInv, m_allNodesFat16B[bvhOffset + leftNodeOffset], trBoxMin, trBoxSize, tNear, pHit->t);
      
      // traversal decision
      //
      leftNodeOffset   = (leftOffsets.z & HIT_LEFT)   != 0 ? leftOffsets.x : leftOffsets.y; // GPU style branch
      if (hitBoth(leftOffsets.z))
      {
        leftNodeOffset = (leftOffsets.z & LEFT_FIRST) != 0 ? leftOffsets.x : leftOffsets.y; // GPU style branch
        stack[top]     = (leftOffsets.z & LEFT_FIRST) != 0 ? leftOffsets.y : leftOffsets.x; // GPU style branch
        top++;
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
      }

      if (missBoth(leftOffsets.z)) // both miss, stack.pop()
      {
        top--;
        leftNodeOffset = stack[std::max(top,0)];
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
      }

    } // end while (searchingForLeaf)

    // leaf node, intersect triangles
    //
    if (top >= 0 && leftNodeOffset != 0xFFFFFFFF)
    {
      const uint32_t start = EXTRACT_START(leftNodeOffset);
      const uint32_t count = EXTRACT_COUNT(leftNodeOffset);

      IntersectAllPrimitivesInLeaf(ray_pos, ray_dir, tNear, instId, geomId, start, count,
                                   pHit);

      #ifdef ENABLE_METRICS
      m_stats.LC++;
      m_stats.LC2++;
      m_stats.TC+=count;
      m_stats.BLB+=count*(3*sizeof(uint32_t) + 3*sizeof(float4));
      #endif
    }

    // continue BVH traversal
    //
    top--;
    leftNodeOffset = stack[std::max(top,0)];
  
    #ifdef ENABLE_METRICS
    m_stats.SOC++;
    m_stats.SBL+=sizeof(uint32_t); 
    #endif
  } // end while (top >= 0)

}

CRT_Hit BVH2FatRTCompact::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
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

bool BVH2FatRTCompact::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  dirAndFar.w *= -1.0f;
  CRT_Hit hit = RayQuery_NearestHit(posAndNear, dirAndFar);
  return (hit.geomId != uint32_t(-1));
}

