#include <iostream>
#include <queue>
#include <stack>
#include <limits>
#include <cassert>

#include <algorithm>
#include <unordered_map>

#include "cbvh_core.h"
#include "cbvh_internal.h"

using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;
using LiteMath::float4x4;

using LiteMath::int2;
using LiteMath::uint2;
using LiteMath::int4;
using LiteMath::uint4;

using LiteMath::min;
using LiteMath::max;
using LiteMath::to_uint32;

using LiteMath::length3f;
using LiteMath::dot3f;
using LiteMath::cross;
using LiteMath::uint;

using cbvh_internal::Box4f;
using cbvh_internal::Triangle4f;

cbvh_internal::RayTriangleHit cbvh_internal::IntersectAllPrimitivesInLeaf(const float4 rayPosAndNear, const float4 rayDirAndFar,
                                                                          const uint* a_indices, uint a_start, uint a_count, const float4* a_vert)
{
  //cvex4::set_ftz();
  const float tNear    = rayPosAndNear[3];

  const float4& ray_pos = rayPosAndNear;
  const float4& ray_dir = rayDirAndFar;

  RayTriangleHit result;
  result.t      = rayDirAndFar[3]; 
  result.primId = -1;
  
  const uint triAddressEnd = a_start + a_count;
 
  for (uint triAddress = a_start; triAddress < triAddressEnd; triAddress += 3)
  { 
    const uint A = a_indices[triAddress + 0];
    const uint B = a_indices[triAddress + 1];
    const uint C = a_indices[triAddress + 2];

    const float4 A_pos = a_vert[A];
    const float4 B_pos = a_vert[B];
    const float4 C_pos = a_vert[C];

    const float4 edge1 = B_pos - A_pos;
    const float4 edge2 = C_pos - A_pos;
    const float4 pvec  = cross(ray_dir, edge2);
    const float4 tvec  = ray_pos - A_pos;
    const float4 qvec  = cross(tvec, edge1);
    const float invDet = 1.0f / std::max(dot3f(edge1, pvec), 1e-6f);

    const float v = dot3f(tvec, pvec)*invDet;
    const float u = dot3f(qvec, ray_dir)*invDet;
    const float t = dot3f(edge2, qvec)*invDet;

    if (v > -1e-6f && u > -1e-6f && (u + v < 1.0f + 1e-6f) && t > tNear && t < result.t)
    {
      result.t      = t;
      result.primId = triAddress/3;
    }
  }

  return result;
}

float4 SafeInverse(float4 d)
{
  const float ooeps = 1.0e-36f; // Avoid div by zero.

  float4 res;
  res.x = 1.0f / (fabs(d.x) > ooeps ? d.x : copysign(ooeps, d.x));
  res.y = 1.0f / (fabs(d.y) > ooeps ? d.y : copysign(ooeps, d.y));
  res.z = 1.0f / (fabs(d.z) > ooeps ? d.z : copysign(ooeps, d.z));
  res.w = d.w;
  return res;
}

#ifndef MAXFLOAT
  #define MAXFLOAT 1e37f
#endif

#define STACK_SIZE 80

using cbvh::BVHNode;
using cbvh_internal::Box4f;

struct NodeSortData
{
  uint leftOffset;
  uint leftOffsetOld;
  float hitDist;
};

struct SortByDist
{
  inline bool operator()(const NodeSortData& struct1, const NodeSortData& struct2)
  {
    return (struct1.hitDist < struct2.hitDist);
  }
};

cbvh_internal::TraversalRes cbvh_internal::BVH4RayTraversalExample(float4 rayPosAndNear, float4 rayDirAndFar,
                                                                   const cbvh::BVHTree& a_tree, const float4* a_vertices)
{
  //cvex4::set_ftz();
  const float4 invDirAndFar = SafeInverse(rayDirAndFar);
  const float tNear = rayPosAndNear[3];

  cbvh_internal::TraversalRes result;
  result.t      = rayDirAndFar[3]; 
  result.primId = -1;

  const Box4f* nodes    = (const Box4f*)a_tree.nodes.data();
  const auto* intervals = a_tree.intervals.data();
  const uint* indices   = a_tree.indicesReordered.data();

  int2  stackData[STACK_SIZE];
  int2* stack  = stackData+2;

  int top                = 0;
  uint leftNodeOffset    = 1;
  uint leftNodeOffsetOld = leftNodeOffset;
  bool searchingForLeaf  = true;

  while (top >= 0)
  {
    while (searchingForLeaf)
    {
      result.numBoxTests += 4;

      const auto node0 = nodes[leftNodeOffset + 0]; 
      const auto node1 = nodes[leftNodeOffset + 1]; 
      const auto node2 = nodes[leftNodeOffset + 2]; 
      const auto node3 = nodes[leftNodeOffset + 3]; 

      const float2 tm0 = LiteMath::Ray4fBox4fIntersection(rayPosAndNear, invDirAndFar, node0.boxMin, node0.boxMax);
      const float2 tm1 = LiteMath::Ray4fBox4fIntersection(rayPosAndNear, invDirAndFar, node1.boxMin, node1.boxMax);
      const float2 tm2 = LiteMath::Ray4fBox4fIntersection(rayPosAndNear, invDirAndFar, node2.boxMin, node2.boxMax);
      const float2 tm3 = LiteMath::Ray4fBox4fIntersection(rayPosAndNear, invDirAndFar, node3.boxMin, node3.boxMax);

      NodeSortData children[4];

      children[0].leftOffset = extractIntW(node0.boxMin);
      children[1].leftOffset = extractIntW(node1.boxMin);
      children[2].leftOffset = extractIntW(node2.boxMin);
      children[3].leftOffset = extractIntW(node3.boxMin);

      children[0].leftOffsetOld = leftNodeOffset + 0;
      children[1].leftOffsetOld = leftNodeOffset + 1;
      children[2].leftOffsetOld = leftNodeOffset + 2;
      children[3].leftOffsetOld = leftNodeOffset + 3;
           
      const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= tNear) && (tm0.x <= result.t);
      const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= tNear) && (tm1.x <= result.t);
      const bool hitChild2 = (tm2.x <= tm2.y) && (tm2.y >= tNear) && (tm2.x <= result.t);
      const bool hitChild3 = (tm3.x <= tm3.y) && (tm3.y >= tNear) && (tm3.x <= result.t);

      children[0].hitDist = hitChild0 ? tm0.x : MAXFLOAT;
      children[1].hitDist = hitChild1 ? tm1.x : MAXFLOAT;
      children[2].hitDist = hitChild2 ? tm2.x : MAXFLOAT;
      children[3].hitDist = hitChild3 ? tm3.x : MAXFLOAT;

      std::sort(children, children+4, SortByDist());

      const bool stackHaveSpace = (top < STACK_SIZE);

      // push/pop stack games
      //
      if (children[3].hitDist < MAXFLOAT && stackHaveSpace)
      {
        stack[top].x = children[3].leftOffset;
        stack[top].y = children[3].leftOffsetOld; 
        top++;
      }

      if (children[2].hitDist < MAXFLOAT && stackHaveSpace)
      {
        stack[top].x = children[2].leftOffset;
        stack[top].y = children[2].leftOffsetOld;
        top++;
      }

      if (children[1].hitDist < MAXFLOAT && stackHaveSpace)
      {
        stack[top].x = children[1].leftOffset;
        stack[top].y = children[1].leftOffsetOld;
        top++;
      }
      
      if (children[0].hitDist < MAXFLOAT)
      {
        leftNodeOffset    = children[0].leftOffset;
        leftNodeOffsetOld = children[0].leftOffsetOld;
      }
      else if (top >= 0)
      {
        top--;
        leftNodeOffset    = stack[top].x;
        leftNodeOffsetOld = stack[top].y;    
      }

      //searchingForLeaf = (leftNodeOffset != 0xFFFFFFFF) && (top >= 0);
      searchingForLeaf = (int(leftNodeOffset) >= 0) && (top >= 0);

    } //  while (searchingForLeaf)

   
    // leaf node, intersect triangles
    //
    if (top >= 0)
    {
      const auto startCount = intervals[leftNodeOffsetOld];
      const auto localHit   = cbvh_internal::IntersectAllPrimitivesInLeaf(rayPosAndNear, rayDirAndFar, 
                                                                          indices, startCount.start*3, startCount.count*3, a_vertices);
      
      if(localHit.t > tNear && localHit.t < result.t)
      {
        result.t       = localHit.t;
        result.primId  = localHit.primId;
        rayDirAndFar.w = localHit.t;
      }      

      result.numPrimTests += startCount.count;   
      result.numLeafesTests++;              
    }

    // pop next node from stack
    //
    top--;                                                  
    leftNodeOffset    = stack[top].x;                         
    leftNodeOffsetOld = stack[top].y;                     
    //searchingForLeaf  = (leftNodeOffset != 0xFFFFFFFF);     
    searchingForLeaf = (int(leftNodeOffset) >= 0);

  } // while (top >= 0)

  return result;
}

cbvh_internal::TraversalRes cbvh_internal::BVH4RayTraversalExample(float4 rayPosAndNear, float4 rayDirAndFar,
                                                                   const cbvh::BVHTree& a_tree, const float4* a_vertices, std::vector<uint32_t>& a_leaves)
{
  //cvex4::set_ftz();
  const float4 invDirAndFar = SafeInverse(rayDirAndFar);
  const float tNear = rayPosAndNear[3];

  cbvh_internal::TraversalRes result;
  result.t      = rayDirAndFar[3]; 
  result.primId = -1;

  const Box4f* nodes    = (const Box4f*)a_tree.nodes.data();
  const auto* intervals = a_tree.intervals.data();
  const uint* indices   = a_tree.indicesReordered.data();

  int2  stackData[STACK_SIZE];
  int2* stack  = stackData+2;

  int top                = 0;
  uint leftNodeOffset    = 1;
  uint leftNodeOffsetOld = leftNodeOffset;
  bool searchingForLeaf  = true;

  a_leaves.reserve(STACK_SIZE);

  while (top >= 0)
  {
    while (searchingForLeaf)
    {
      result.numBoxTests += 4;

      const auto node0 = nodes[leftNodeOffset + 0]; 
      const auto node1 = nodes[leftNodeOffset + 1]; 
      const auto node2 = nodes[leftNodeOffset + 2]; 
      const auto node3 = nodes[leftNodeOffset + 3]; 

      const float2 tm0 = LiteMath::Ray4fBox4fIntersection(rayPosAndNear, invDirAndFar, node0.boxMin, node0.boxMax);
      const float2 tm1 = LiteMath::Ray4fBox4fIntersection(rayPosAndNear, invDirAndFar, node1.boxMin, node1.boxMax);
      const float2 tm2 = LiteMath::Ray4fBox4fIntersection(rayPosAndNear, invDirAndFar, node2.boxMin, node2.boxMax);
      const float2 tm3 = LiteMath::Ray4fBox4fIntersection(rayPosAndNear, invDirAndFar, node3.boxMin, node3.boxMax);

      NodeSortData children[4];

      children[0].leftOffset = extractIntW(node0.boxMin);
      children[1].leftOffset = extractIntW(node1.boxMin);
      children[2].leftOffset = extractIntW(node2.boxMin);
      children[3].leftOffset = extractIntW(node3.boxMin);

      children[0].leftOffsetOld = leftNodeOffset + 0;
      children[1].leftOffsetOld = leftNodeOffset + 1;
      children[2].leftOffsetOld = leftNodeOffset + 2;
      children[3].leftOffsetOld = leftNodeOffset + 3;
           
      const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= tNear) && (tm0.x <= result.t);
      const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= tNear) && (tm1.x <= result.t);
      const bool hitChild2 = (tm2.x <= tm2.y) && (tm2.y >= tNear) && (tm2.x <= result.t);
      const bool hitChild3 = (tm3.x <= tm3.y) && (tm3.y >= tNear) && (tm3.x <= result.t);

      children[0].hitDist = hitChild0 ? tm0.x : MAXFLOAT;
      children[1].hitDist = hitChild1 ? tm1.x : MAXFLOAT;
      children[2].hitDist = hitChild2 ? tm2.x : MAXFLOAT;
      children[3].hitDist = hitChild3 ? tm3.x : MAXFLOAT;

      std::sort(children, children+4, SortByDist());

      const bool stackHaveSpace = (top < STACK_SIZE);

      // push/pop stack games
      //
      if (children[3].hitDist < MAXFLOAT && stackHaveSpace)
      {
        stack[top].x = children[3].leftOffset;
        stack[top].y = children[3].leftOffsetOld; 
        top++;
      }

      if (children[2].hitDist < MAXFLOAT && stackHaveSpace)
      {
        stack[top].x = children[2].leftOffset;
        stack[top].y = children[2].leftOffsetOld;
        top++;
      }

      if (children[1].hitDist < MAXFLOAT && stackHaveSpace)
      {
        stack[top].x = children[1].leftOffset;
        stack[top].y = children[1].leftOffsetOld;
        top++;
      }
      
      if (children[0].hitDist < MAXFLOAT)
      {
        leftNodeOffset    = children[0].leftOffset;
        leftNodeOffsetOld = children[0].leftOffsetOld;
      }
      else if (top >= 0)
      {
        top--;
        leftNodeOffset    = stack[top].x;
        leftNodeOffsetOld = stack[top].y;    
      }

      //searchingForLeaf = (leftNodeOffset != 0xFFFFFFFF) && (top >= 0);
      searchingForLeaf = (int(leftNodeOffset) >= 0) && (top >= 0);

    } //  while (searchingForLeaf)

   
    // leaf node, intersect triangles
    //
    if (top >= 0)
    {
      a_leaves.push_back(leftNodeOffsetOld);
      const auto startCount = intervals[leftNodeOffsetOld];
      const auto localHit   = cbvh_internal::IntersectAllPrimitivesInLeaf(rayPosAndNear, rayDirAndFar, 
                                                                          indices, startCount.start*3, startCount.count*3, a_vertices);
      
      if(localHit.t > tNear && localHit.t < result.t)
      {
        result.t       = localHit.t;
        result.primId  = localHit.primId;
        rayDirAndFar.w = localHit.t;
      }      

      result.numPrimTests += startCount.count;   
      result.numLeafesTests++;              
    }

    // pop next node from stack
    //
    top--;                                                  
    leftNodeOffset    = stack[top].x;                         
    leftNodeOffsetOld = stack[top].y;                     
    //searchingForLeaf  = (leftNodeOffset != 0xFFFFFFFF);     
    searchingForLeaf = (int(leftNodeOffset) >= 0);

  } // while (top >= 0)

  return result;
}
