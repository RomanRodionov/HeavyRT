#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>

#include "BVH4HalfRT.h"

void BVH4HalfRT::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir, 
                                              float tNear, uint32_t instId, uint32_t geomId,
                                              uint32_t a_start, uint32_t a_count,
                                              CRT_Hit* pHit)
{
  const uint2 a_geomOffsets = m_geomOffsets[geomId]; //uint2(m_geomOffsets[geomId].x,0); //m_geomOffsets[geomId];

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
    const float3 pvec  = cross(ray_dir, edge2);
    const float3 tvec  = ray_pos - A_pos;
    const float3 qvec  = cross(tvec, edge1);
    
    const float invDet = 1.0f / dot(edge1, pvec);
    const float v      = dot(tvec, pvec)*invDet;
    const float u      = dot(qvec, ray_dir)*invDet;
    const float t      = dot(edge2, qvec)*invDet;
  
    if (v >= -1e-6f && u >= -1e-6f && (u + v <= 1.0f+1e-6f) && t > tNear && t < pHit->t) // if (v > -1e-6f && u > -1e-6f && (u + v < 1.0f+1e-6f) && t > tMin && t < hit.t)
    {
      pHit->t         = t;
      pHit->primId    = triId;
      pHit->instId    = instId;
      pHit->geomId    = geomId;  
      pHit->coords[0] = u;
      pHit->coords[1] = v;
    }
  }
}

static inline float2 unpackFloat2x16(uint32_t v)
{
  LiteMath::half xy[2];
  memcpy(xy, &v, sizeof(uint32_t));
  return float2(xy[0], xy[1]);
}

static inline BVHNode UnpackBVHNode(uint4 data)
{
  const float2 XY  = float2(unpackFloat2x16(data.x));
  const float2 ZX  = float2(unpackFloat2x16(data.y));
  const float2 YZ  = float2(unpackFloat2x16(data.z));
  BVHNode node;
  node.boxMin.x    = XY.x;
  node.boxMin.y    = XY.y;
  node.boxMin.z    = ZX.x;
  node.boxMax.x    = ZX.y;
  node.boxMax.y    = YZ.x;
  node.boxMax.z    = YZ.y;
  node.leftOffset  = data.w;
  node.escapeIndex = 0;
  return node;
}

void BVH4HalfRT::BVH4Traverse(const float3 ray_pos, const float3 ray_dir, float tNear, uint32_t instId, uint32_t geomId, uint32_t stack[STACK_SIZE],
                              bool stopOnFirstHit, CRT_Hit* pHit)

{
  const uint32_t bvhOffset = m_bvhOffsets[geomId];
  uint32_t leftNodeOffset  = 0; 
  const float3 rayDirInv   = SafeInverse(ray_dir);
  int top                  = 0;

  while (top >= 0 && !(stopOnFirstHit && pHit->primId != uint32_t(-1)))
  {
    bool needStackPop = false;
    if ((leftNodeOffset & LEAF_BIT) == 0)
    {
      #ifdef ENABLE_METRICS
      m_stats.NC  += 4;
      m_stats.BLB += 4*sizeof(int);
      #endif
      
      const BVHNode node0 = UnpackBVHNode(m_allNodes16[bvhOffset + leftNodeOffset + 0]); 
      const BVHNode node1 = UnpackBVHNode(m_allNodes16[bvhOffset + leftNodeOffset + 1]); 
      const BVHNode node2 = UnpackBVHNode(m_allNodes16[bvhOffset + leftNodeOffset + 2]); 
      const BVHNode node3 = UnpackBVHNode(m_allNodes16[bvhOffset + leftNodeOffset + 3]); 

      const float2 tm0 = RayBoxIntersection2(ray_pos, rayDirInv, node0.boxMin, node0.boxMax);
      const float2 tm1 = RayBoxIntersection2(ray_pos, rayDirInv, node1.boxMin, node1.boxMax);
      const float2 tm2 = RayBoxIntersection2(ray_pos, rayDirInv, node2.boxMin, node2.boxMax);
      const float2 tm3 = RayBoxIntersection2(ray_pos, rayDirInv, node3.boxMin, node3.boxMax);

      const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= tNear) && (tm0.x <= pHit->t);
      const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= tNear) && (tm1.x <= pHit->t);
      const bool hitChild2 = (tm2.x <= tm2.y) && (tm2.y >= tNear) && (tm2.x <= pHit->t);
      const bool hitChild3 = (tm3.x <= tm3.y) && (tm3.y >= tNear) && (tm3.x <= pHit->t);
      
      // todo: try to join int and float to single float2 or int2 element to make less comparison operators!
      //float2 n0 = float2(hitChild0 ? tm0.x : MAXFLOAT, as_float(node0.leftOffset));
      //float2 n1 = float2(hitChild1 ? tm1.x : MAXFLOAT, as_float(node1.leftOffset));
      //float2 n2 = float2(hitChild2 ? tm2.x : MAXFLOAT, as_float(node2.leftOffset));
      //float2 n3 = float2(hitChild3 ? tm3.x : MAXFLOAT, as_float(node3.leftOffset));
    
      int4 children  = int4(node0.leftOffset, node1.leftOffset, node2.leftOffset, node3.leftOffset);
      float4 hitMinD = float4(hitChild0 ? tm0.x : MAXFLOAT,
                              hitChild1 ? tm1.x : MAXFLOAT,
                              hitChild2 ? tm2.x : MAXFLOAT,
                              hitChild3 ? tm3.x : MAXFLOAT);

      // sort tHit and children
      //
      const bool lessXY = (hitMinD.y < hitMinD.x);
      const bool lessWZ = (hitMinD.w < hitMinD.z);
      {
        const int   w_childrenX = lessXY ? children.y : children.x;
        const int   w_childrenY = lessXY ? children.x : children.y;
        const float w_hitTimesX = lessXY ? hitMinD.y  : hitMinD.x;
        const float w_hitTimesY = lessXY ? hitMinD.x  : hitMinD.y;

        const int   w_childrenZ = lessWZ ? children.w : children.z;
        const int   w_childrenW = lessWZ ? children.z : children.w;
        const float w_hitTimesZ = lessWZ ? hitMinD.w  : hitMinD.z;
        const float w_hitTimesW = lessWZ ? hitMinD.z  : hitMinD.w;

        children.x = w_childrenX;
        children.y = w_childrenY;
        hitMinD.x  = w_hitTimesX;
        hitMinD.y  = w_hitTimesY;

        children.z = w_childrenZ;
        children.w = w_childrenW;
        hitMinD.z  = w_hitTimesZ;
        hitMinD.w  = w_hitTimesW;
      }

      const bool lessZX = (hitMinD.z < hitMinD.x);
      const bool lessWY = (hitMinD.w < hitMinD.y);
      {
        const int   w_childrenX = lessZX ? children.z : children.x;
        const int   w_childrenZ = lessZX ? children.x : children.z;
        const float w_hitTimesX = lessZX ? hitMinD.z  : hitMinD.x;
        const float w_hitTimesZ = lessZX ? hitMinD.x  : hitMinD.z;

        const int   w_childrenY = lessWY ? children.w : children.y;
        const int   w_childrenW = lessWY ? children.y : children.w;
        const float w_hitTimesY = lessWY ? hitMinD.w  : hitMinD.y;
        const float w_hitTimesW = lessWY ? hitMinD.y  : hitMinD.w;

        children.x = w_childrenX;
        children.z = w_childrenZ;
        hitMinD.x  = w_hitTimesX; 
        hitMinD.z  = w_hitTimesZ; 

        children.y = w_childrenY;
        children.w = w_childrenW;
        hitMinD.y  = w_hitTimesY;
        hitMinD.w  = w_hitTimesW;
      }

      const bool lessZY = (hitMinD.z < hitMinD.y);
      {
        const int   w_childrenY = lessZY ? children.z : children.y;
        const int   w_childrenZ = lessZY ? children.y : children.z;
        const float w_hitTimesY = lessZY ? hitMinD.z  : hitMinD.y;
        const float w_hitTimesZ = lessZY ? hitMinD.y  : hitMinD.z;

        children.y = w_childrenY;
        children.z = w_childrenZ;
        hitMinD.y  = w_hitTimesY;
        hitMinD.z  = w_hitTimesZ;
      }

      // push/pop stack games
      //
      if (hitMinD.w < MAXFLOAT)
      {
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
        stack[top] = children.w;
        top++;
      }

      if (hitMinD.z < MAXFLOAT)
      {
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
        stack[top] = children.z;
        top++;
      }

      if (hitMinD.y < MAXFLOAT)
      {
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
        stack[top] = children.y;
        top++;
      }

      if (hitMinD.x < MAXFLOAT)
      {
        #ifdef ENABLE_METRICS
        m_stats.SOC++;
        m_stats.SBL+=sizeof(uint32_t); 
        #endif
        leftNodeOffset = children.x;
      }
      needStackPop = (hitMinD.x >= MAXFLOAT);
                        
    } // end if(notLeaf)
   
    // leaf node, intersect triangles
    //
    if ((leftNodeOffset & LEAF_BIT) != 0 && leftNodeOffset != 0xFFFFFFFF) 
    {
      const uint32_t start = EXTRACT_START(leftNodeOffset);
      const uint32_t count = EXTRACT_COUNT(leftNodeOffset);
      IntersectAllPrimitivesInLeaf(ray_pos, ray_dir, tNear, instId, geomId, start, count, pHit);
      needStackPop = true;
      #ifdef ENABLE_METRICS
      m_stats.LC++;
      m_stats.LC2++;
      m_stats.TC  += count;
      m_stats.BLB += count*(3*sizeof(uint32_t) + 3*sizeof(float4));
      #endif          
    }

    // pop next node from stack
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

CRT_Hit BVH4HalfRT::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  bool stopOnFirstHit = (dirAndFar.w <= 0.0f);
  if(stopOnFirstHit)
    dirAndFar.w *= -1.0f;

  #ifdef ENABLE_METRICS
  m_stats.raysNumber++;
  ResetVarLC();
  #endif

  CRT_Hit hit;
  hit.t      = dirAndFar.w;
  hit.primId = uint32_t(-1);
  hit.instId = uint32_t(-1);
  hit.geomId = uint32_t(-1);

  const float3 rayDirInv = SafeInverse(to_float3(dirAndFar));

  uint32_t stack[STACK_SIZE];
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
      m_stats.BLB += 1*sizeof(BVHNode);
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
  
      BVH4Traverse(ray_pos, ray_dir, posAndNear.w, instId, geomId, stack, stopOnFirstHit, &hit);
    }

  } while (nodeIdx < 0xFFFFFFFE && !(stopOnFirstHit && hit.primId != uint32_t(-1)));
  

  #ifdef REMAP_PRIM_ID
  if(hit.geomId < uint32_t(-1)) 
  {
    const uint2 geomOffsets = m_geomOffsets[hit.geomId];
    hit.primId = m_primIndices[geomOffsets.x/3 + hit.primId];
  }
  #endif
 
  return hit;
}

bool BVH4HalfRT::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  dirAndFar.w *= -1.0f;
  CRT_Hit hit = RayQuery_NearestHit(posAndNear, dirAndFar);
  return (hit.geomId != uint32_t(-1));
}
