#include <iostream>
#include <cassert>

#include "cbvh_core.h"
#include "cbvh_internal.h"
#include "bvhree_tests.h"

#include <cstdint>
#include <cstdlib> // int rand();

float inline rnd()
{
  return (rand() % RAND_MAX) / (float)RAND_MAX;
}

using LiteMath::float4;
using LiteMath::Box4f;

void RandomRaysInBox(const LiteMath::Box4f& scnBox, int a_raysNum,
                     LiteMath::float4* raysPos, LiteMath::float4* raysDir)
{
  const float4 boxSize = scnBox.boxMax - scnBox.boxMin;
  
  for(int i=0;i<a_raysNum;i++)
  {
    const float4 randOffset = 2.0f*LiteMath::normalize(float4(rnd(),rnd(),rnd(),0.0f)) - float4(1,1,1,0); // [0,1] --> [-1,1]
    raysPos[i]   = scnBox.boxMin + float4(rnd(),rnd(),rnd(),0.0f)*boxSize; // #TODO: use random sphere direction
    raysDir[i]   = LiteMath::normalize(randOffset);                        // #TODO: use random sphere direction 
    raysPos[i].w = 0.0f;
    raysDir[i].w = LiteMath::INF_POSITIVE;
  }
}


void RandomRaysOutOfBox(const LiteMath::Box4f& scnBox, int a_raysNum,
                        LiteMath::float4* raysPos, LiteMath::float4* raysDir)
{
  const float4 boxSize     = scnBox.boxMax - scnBox.boxMin;
  const float sphereRadius = 0.5f*LiteMath::length3f(boxSize);
  const float4 boxCenter   = 0.5f*(scnBox.boxMax + scnBox.boxMin);

  const float4 randOffset    = 2.0f*LiteMath::normalize(float4(rnd(),rnd(),rnd(),0.0f)) - float4(1,1,1,0); // [0,1] --> [-1,1]
  const float4 pointInTheBox = scnBox.boxMin + float4(rnd(),rnd(),rnd(),0.0f)*boxSize;

  for(int i=0;i<a_raysNum;i++)
  {
    raysPos[i]   = boxCenter + randOffset*2.0f*sphereRadius;
    raysDir[i]   = LiteMath::normalize(pointInTheBox - raysPos[i]);
    raysPos[i].w = 0.0f;
    raysDir[i].w = 4.0f*sphereRadius;
  }
}

void cbvh_test::GenerateUniformRandomRays(const LiteMath::Box4f& scnBox, int TST_RAYS,
                                          LiteMath::float4* raysPos, LiteMath::float4* raysDir)
{
  RandomRaysInBox   (scnBox, (TST_RAYS/2), raysPos, raysDir);
  RandomRaysOutOfBox(scnBox, (TST_RAYS/2), raysPos+(TST_RAYS/2), raysDir+(TST_RAYS/2));
}


void GenerateTestRandomRays(const Box4f& scnBox, std::vector<float4>& raysPos, std::vector<float4>& raysDir, int TST_RAYS = 512)
{
  raysPos.resize(TST_RAYS);
  raysDir.resize(TST_RAYS);

  cbvh_test::GenerateUniformRandomRays(scnBox, TST_RAYS,
                                       raysPos.data(), raysDir.data());
}

bool cbvh_test::TestRayTraversal(const cbvh::BVHTree& in_bvhTree, 
                                 const LiteMath::float4* in_vPos4f, size_t a_maxVertexCount,
                                 const uint* a_indicesOriginal, size_t a_indexNum)
{
  if(in_bvhTree.nodes.size() == 0)
    return false;

  const Box4f* boxData = (const Box4f*)in_bvhTree.nodes.data();

  Box4f scnBox;
  for(int i=0;i<a_indexNum;i+=3)
  {
    const auto A = a_indicesOriginal[i+0];
    const auto B = a_indicesOriginal[i+1];
    const auto C = a_indicesOriginal[i+2]; 

    assert(A < a_maxVertexCount);
    assert(B < a_maxVertexCount);
    assert(C < a_maxVertexCount);

    scnBox.include(in_vPos4f[A]);
    scnBox.include(in_vPos4f[B]);
    scnBox.include(in_vPos4f[C]);
  }

  Box4f bvhBox = boxData[0];

  const float diff1 = length3f(scnBox.boxMin - bvhBox.boxMin);
  const float diff2 = length3f(scnBox.boxMax - bvhBox.boxMax);

  if(diff1 > 1e-6f || diff2 > 1e-6f)
    return false;

  std::vector<float4> raysPos, raysDir;
  //GenerateTestRandomRays(scnBox, raysPos, raysDir, 256);
  raysPos.emplace_back(-0.00264156517, -0.00264156517, 0.00543617224, 0);
  raysDir.emplace_back(-0.203098252, 0.733451784, -0.648690641, +INFINITY);

  const float4 boxSize   = scnBox.boxMax - scnBox.boxMin;
  const float4 boxCenter = 0.5f*(scnBox.boxMax + scnBox.boxMin);
  const float sphereRadius = 0.5f*LiteMath::length3f(boxSize);

  bool failed = false;
  for(int i=0;i<raysPos.size();i++)
  {
    auto hitBvh = cbvh_internal::BVH4RayTraversalExample(raysPos[i], raysDir[i], in_bvhTree, in_vPos4f); 
    auto hitBF  = cbvh_internal::IntersectAllPrimitivesInLeaf(raysPos[i], raysDir[i], a_indicesOriginal, 0, a_indexNum, in_vPos4f);
    //auto hitBF  = cbvh_internal::IntersectAllPrimitivesInLeaf(raysPos[i], raysDir[i], in_bvhTree.indicesReordered.data(), 0, in_bvhTree.indicesReordered.size(), in_vPos4f);

    if(fabs(hitBvh.t - hitBF.t) > 1e-5f)
    {
      failed = true;
      std::cout << std::endl;
      std::cout << "[TestRayTraversal]: i = " << i << "; t1 = " << hitBvh.t << "; t2 = " << hitBF.t << std::endl;
      for(int j=0;j<4;j++)
        std::cout << raysPos[i][j] << " ";
      std::cout << std::endl;
      for(int j=0;j<4;j++)
        std::cout << raysDir[i][j] << " ";
      std::cout << std::endl;
      std::cout.flush();
      break;
    }

    if(i%16 == 0)
    {
      std::cout << "[TestRayTraversal]: progress = " << 100.0f*(float(i)/float(raysPos.size())) << "%  \r";
      std::cout.flush();
    }
  }
  std::cout << std::endl;

  return !failed;
}

cbvh_test::TraversalMetricsData cbvh_test::ComputeRayTraversalMetrics(const cbvh::BVHTree& in_bvhTree, 
                                                                      const LiteMath::float4* in_vPos4f, size_t a_maxVertexCount,
                                                                      const uint* a_indicesOriginal, size_t a_indexNum)
{
  if(in_bvhTree.nodes.size() == 0)
    return TraversalMetricsData();

  const Box4f* boxData = (const Box4f*)in_bvhTree.nodes.data();
  const Box4f scnBox   = boxData[0];

  std::vector<float4> raysPos, raysDir;
  GenerateTestRandomRays(scnBox, raysPos, raysDir, 16384);

  int64_t totalNodes = 0;
  int64_t totalLeafs = 0;
  int64_t totalPrims = 0;
  double  totalLeafesSquare = 0.0;

  for(int i=0;i<raysPos.size();i++)
  {
    auto hitBvh = cbvh_internal::BVH4RayTraversalExample(raysPos[i], raysDir[i], in_bvhTree, in_vPos4f); 
    totalNodes += hitBvh.numBoxTests;
    totalLeafs += hitBvh.numLeafesTests;
    totalPrims += hitBvh.numPrimTests;
    totalLeafesSquare += double(hitBvh.numLeafesTests*hitBvh.numLeafesTests);
  }

  float avgLC       = float(totalLeafs)/float(raysPos.size());
  float avgLCSquare = float(totalLeafesSquare)/float(raysPos.size());

  assert(avgLCSquare - avgLC*avgLC > 0.0f);

  TraversalMetricsData res;
  res.avgNodesCount  = float(totalNodes)/float(raysPos.size());
  res.avgLeafesCount = float(totalLeafs)/float(raysPos.size());
  res.avgPrimsCount  = float(totalPrims)/float(raysPos.size());
  res.varLeafesCount = sqrtf(avgLCSquare - avgLC*avgLC);
  return res;
}

std::vector<Box4f> cbvh_test::TriangleBoxesFromMesh(const float4* a_vpos, size_t a_vertNum, const uint* a_indices, size_t a_indexNum)
{
  std::vector<Box4f> boxes(a_indexNum/3);

  const int trisNum = int(a_indexNum/3);
  for(int i=0;i<trisNum;i++)
  {
    for(int j=0;j<3;j++)
    {
      const float4 vert = a_vpos[a_indices[i*3+j]];
      boxes[i].boxMin = LiteMath::min(boxes[i].boxMin, vert);
      boxes[i].boxMax = LiteMath::max(boxes[i].boxMax, vert);
    }

    boxes[i].setStart(i);
    boxes[i].setCount(1);
  }

  return boxes;
}
