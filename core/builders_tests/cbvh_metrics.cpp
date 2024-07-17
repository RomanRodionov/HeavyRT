#include <iostream>
#include <queue>
#include <stack>
#include <algorithm>

#include "cbvh_metrics.h"
#include "cbvh_internal.h"

float cbvh_metrics::MetricSAH(const BVHTree &in_bvhTree, const float4* in_verts, unsigned int in_verts_count,
                              const float in_innerC, const float in_leafC) {
  float SAHsum = 0.0;
  // iterate over all nodes
  auto nodeDiv = 1.0f/ in_bvhTree.nodes.size();
  for(int i = 0; i < in_bvhTree.nodes.size(); ++i){
    // calculate SAH with each child node
    auto& currentNode = in_bvhTree.nodes[i];
    if(currentNode.leftOffset != -3 && currentNode.leftOffset != -1){
      for(int child = 0; child < 4; ++child){
        auto& currentChild = in_bvhTree.nodes[currentNode.leftOffset + child];
        float currentCost = 0.0f;
        // in case of the leaf
        if (currentChild.leftOffset == -1)
          currentCost = in_leafC * in_bvhTree.intervals[currentNode.leftOffset + child].count;
        // in case of empty node
        else
        if (currentChild.leftOffset == -3)
          continue;
        // in case of regular node
        else
          currentCost = in_innerC;

        auto surf1 = cbvh_internal::SurfaceArea(currentNode.boxMin,  currentNode.boxMax);
        auto surf2 = cbvh_internal::SurfaceArea(currentChild.boxMin, currentChild.boxMax);
        if(surf1 == 0 || surf2 == 0)
          continue;
        
        if(currentChild.boxMax[0] >= currentChild.boxMin[0] && 
           currentChild.boxMax[1] >= currentChild.boxMin[1] && 
           currentChild.boxMax[2] >= currentChild.boxMin[2])  // do not process invalid node
        {
          auto ret = currentCost * (cbvh_internal::SurfaceArea(currentChild.boxMin, currentChild.boxMax)/
          cbvh_internal::SurfaceArea(currentNode.boxMin,  currentNode.boxMax)) * nodeDiv * 0.25f;
          if(!std::isfinite(ret)){
            int gg = 4;
            std::cout << gg;
          }

          SAHsum += ret;
        }
      }
    }
  }
  
  //std::cout << "SAH: \t" << SAHsum << std::endl;
  return SAHsum; //*float(in_bvhTree.nodes.size());
}

static inline float SurfaceAreaOfTriangle(float4 in_a, float4 in_b, float4 in_c)
{
  auto u = LiteMath::cross3(in_c - in_a, in_b - in_a);
  return LiteMath::length3f(u) / 2.0;
}

// Calculates the vector of node indices that we have to test the intersection with
inline void InitVectorWithIndices(std::vector<unsigned int>& in_vector, unsigned int in_vectorsize, std::vector<unsigned int> in_throwIndices){
  // Initialize the vector with indices from 0 to in_vectorsize
  for(int i = 0; i < in_vectorsize; ++i){
    in_vector.push_back(i);
  }

  // Sorting is needed to correctly delete indices from larger to smaller indices
  std::sort(in_throwIndices.begin(), in_throwIndices.end());

  // deleting elements by index in reverse
  for(int i = in_throwIndices.size() - 1; i >= 0; --i){
    if(i < in_vector.size() && i >= 0)
      in_vector.erase(in_vector.begin() + i);
  }

  //std::cout << "test";
}

// Naive brute-force EPO metric
float cbvh_metrics::MetricEPOTriBF(const BVHTree &in_bvhTree, const float4* in_verts, unsigned int in_verts_count,
                                 const float in_innerC, const float in_leafC) {
  float EPOsum = 0.0;

  std::vector<unsigned int> nodesPassed;   // save indices of parent nodes

  std::stack<std::pair<unsigned int, unsigned int>> traverseQueue; // queue of {index, depth} used for tree traversing
  traverseQueue.push({0, 0});

  while(traverseQueue.size() != 0){
    auto currentIndex = traverseQueue.top();
    auto& currentNode = in_bvhTree.nodes[currentIndex.first];
    auto& currentInterval = in_bvhTree.intervals[currentIndex.first];

    //correcting node level in case we jumped a few levels back instead of one
    while(nodesPassed.size() != currentIndex.second)
      nodesPassed.pop_back();
    nodesPassed.push_back(currentIndex.first);

    traverseQueue.pop();

    if(cbvh::IsLeaf(currentNode)){
      std::vector<unsigned int> nodesToPass;
      nodesToPass.reserve(in_bvhTree.nodes.size());
      InitVectorWithIndices(nodesToPass, in_bvhTree.nodes.size(), nodesPassed);

      auto triDiv = 1.0 / currentInterval.count;
      for(int tri = currentInterval.start * 3; tri < (currentInterval.start + currentInterval.count) * 3; tri += 3) {

        float4 posA = in_verts[in_bvhTree.indicesReordered[tri]];
        float4 posB = in_verts[in_bvhTree.indicesReordered[tri + 1]];
        float4 posC = in_verts[in_bvhTree.indicesReordered[tri + 2]];

        Box4f tempBox;
        tempBox.include(posA);
        tempBox.include(posB);
        tempBox.include(posC);

        auto tempBoxArea = tempBox.surfaceArea();

        auto nodesDiv = 1.0 / nodesToPass.size();
        for (int j = 0; j < nodesToPass.size(); ++j) {

          auto &nodeJ = in_bvhTree.nodes[nodesToPass[j]];
          if (cbvh::IsEmpty(nodeJ))
            continue;
          Box4f nodeBox;
          nodeBox.boxMin.x = nodeJ.boxMin.x;
          nodeBox.boxMin.y = nodeJ.boxMin.y;
          nodeBox.boxMin.z = nodeJ.boxMin.z;
          nodeBox.boxMax.x = nodeJ.boxMax.x;
          nodeBox.boxMax.y = nodeJ.boxMax.y;
          nodeBox.boxMax.z = nodeJ.boxMax.z;


          auto interBox = BoxBoxOverlap(tempBox, nodeBox);
          if (LiteMath::length3f(interBox.boxMax - interBox.boxMin) < 0.0001)
            continue;
          auto interBoxArea = interBox.surfaceArea();

          auto boxBoxRatio = interBoxArea / tempBoxArea;

          /*
          if(!std::isfinite(boxBoxRatio))
            std::cout << "1";
          if(!std::isfinite(SurfaceAreaOfTriangle(posA, posB, posC)))
            std::cout << "2";
          */

          EPOsum += in_leafC * boxBoxRatio * SurfaceAreaOfTriangle(posA, posB, posC) * nodesDiv * triDiv;
        }
      }
      nodesPassed.pop_back();
    }
    else
    if(cbvh::IsEmpty(currentNode)){
      //do nothing
      //std::cout << "empty\n";
    }
    else{
      for(int j = 0; j < 4; ++j) {
        traverseQueue.push({currentNode.leftOffset + j, currentIndex.second + 1});
      }
    }
  }

  //std::cout << "EPO: \t" << EPOsum << std::endl;
  return EPOsum;
}

bool CheckIfContains(unsigned int in_value, std::vector<unsigned int> in_nodesPassed){
  for(int i = 0; i < in_nodesPassed.size(); ++i){
    if(in_nodesPassed[i] == in_value)
      return true;
  }
  return false;
}

float TestIntersections(const BVHTree &in_bvhTree, Box4f in_box, std::vector<unsigned int> in_nodesPassed,
                        float in_boxMult, const float in_innerC, const float in_leafC){

  auto boxArea = in_box.surfaceArea();

  std::vector<unsigned int> nodes;
  nodes.push_back(0);

  float EPOsum = 0.0f;

  while(nodes.size() != 0) {
    auto &currentNode = in_bvhTree.nodes[nodes.back()];

    if(IsEmpty(currentNode)) {
      nodes.pop_back();
      continue;
    }
    else
    if(CheckIfContains(nodes.back(), in_nodesPassed)) {
      nodes.pop_back();
      if (!IsLeaf(currentNode)) {
        for (int i = 0; i < 4; ++i) {
          nodes.push_back(currentNode.leftOffset + i);
        }
      }
      continue;
    }

    Box4f nodeBox;
    nodeBox.boxMin.x = currentNode.boxMin.x;
    nodeBox.boxMin.y = currentNode.boxMin.y;
    nodeBox.boxMin.z = currentNode.boxMin.z;

    nodeBox.boxMax.x = currentNode.boxMax.x;
    nodeBox.boxMax.y = currentNode.boxMax.y;
    nodeBox.boxMax.z = currentNode.boxMax.z;

    auto interBox = BoxBoxOverlap(in_box, nodeBox);
    auto interBoxArea = interBox.surfaceArea();
    if (interBoxArea < 0.0000000001) {
      nodes.pop_back();
      continue;
    }

    auto boxBoxRatio = interBoxArea / boxArea;

    EPOsum += in_leafC * boxBoxRatio * in_boxMult;

    nodes.pop_back();

    if (!IsLeaf(currentNode)) {
      for(int i = 0; i < 4; ++i){
        nodes.push_back(currentNode.leftOffset + i);
      }
    }
  }
  return EPOsum;
}

float cbvh_metrics::MetricEPOTri(const BVHTree &in_bvhTree, const float4* in_verts,
                                 unsigned int in_verts_count, const float in_innerC, const float in_leafC) {
  float EPOsum = 0.0;

  std::vector<unsigned int> nodesPassed;   // save indices of parent nodes

  std::stack<std::pair<unsigned int, unsigned int>> traverseQueue; // queue of {index, depth} used for tree traversing
  traverseQueue.push({0, 0});

  while(traverseQueue.size() != 0){
    auto currentIndex = traverseQueue.top();
    auto& currentNode = in_bvhTree.nodes[currentIndex.first];
    auto& currentInterval = in_bvhTree.intervals[currentIndex.first];

    //correcting node level in case we jumped a few levels back instead of one
    while(nodesPassed.size() != currentIndex.second)
      nodesPassed.pop_back();
    nodesPassed.push_back(currentIndex.first);

    traverseQueue.pop();

    if(cbvh::IsLeaf(currentNode)){

      //auto triDiv = 1.0 / currentInterval.count;
      for(int tri = currentInterval.start * 3; tri < (currentInterval.start + currentInterval.count) * 3; tri += 3){

        float4 posA = in_verts[in_bvhTree.indicesReordered[tri]];
        float4 posB = in_verts[in_bvhTree.indicesReordered[tri + 1]];
        float4 posC = in_verts[in_bvhTree.indicesReordered[tri + 2]];

        Box4f tempBox;
        tempBox.include(posA);
        tempBox.include(posB);
        tempBox.include(posC);

        auto nodesDiv = 1.0 / (in_bvhTree.nodes.size() - nodesPassed.size());
        auto eqMult = SurfaceAreaOfTriangle(posA, posB, posC);

        EPOsum += TestIntersections(in_bvhTree, tempBox, nodesPassed, eqMult, in_innerC, in_leafC);
      }
      nodesPassed.pop_back();
    }
    else
    if(cbvh::IsEmpty(currentNode)){
      //do nothing
      //std::cout << "empty\n";
    }
    else{
      for(int j = 0; j < 4; ++j) {
        traverseQueue.push({currentNode.leftOffset + j, currentIndex.second + 1});
      }
    }
  }

  double sceneSurfaceArea = 0.0;
  for(int i = 0; i < in_bvhTree.indicesReordered.size(); i += 3){
    float4 posA = in_verts[in_bvhTree.indicesReordered[i]];
    float4 posB = in_verts[in_bvhTree.indicesReordered[i + 1]];
    float4 posC = in_verts[in_bvhTree.indicesReordered[i + 2]];

    sceneSurfaceArea += SurfaceAreaOfTriangle(posA, posB, posC);
  }

  //std::cout << "EPO: \t" << EPOsum << std::endl;
  return EPOsum / sceneSurfaceArea;
}

