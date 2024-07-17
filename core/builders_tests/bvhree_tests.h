#pragma once
#include <stack>
#include <set>
#include <iostream>
#include <cassert>

#include "cbvh_core.h"
#include "cbvh_internal.h"

#include "LiteMath.h"

namespace cbvh_test{

    using LiteMath::uint;
    using LiteMath::float4;
    using LiteMath::dot3f;
    using LiteMath::Box4f;

    /*
     * BVH depth test. Traverses the tree and prints minimum and maximum depth.
     */
    static bool TestBVHdepth(const cbvh::BVHTree& in_bvhTree, int* pMinDepth, int* pMaxDepth)
    {
      if(in_bvhTree.nodes.size() == 1)
        return IsLeaf(in_bvhTree.nodes[0]);

      std::stack<std::pair<cbvh::BVHNode, int> > stack;
      stack.push(std::pair<cbvh::BVHNode, int>(in_bvhTree.nodes[0], 0));
      
      int bvhMaxDepth = 0;
      int bvhMinDepth = 10000;
      while(!stack.empty())
      {
        cbvh::BVHNode currentNode = stack.top().first;
        int currentDepth = stack.top().second;
        stack.pop();
        if(currentNode.leftOffset != 0xFFFFFFFF && currentNode.leftOffset != 0xFFFFFFFD){
          for(int i = 0; i < 4; ++i) {
            stack.push(std::pair<cbvh::BVHNode, int>(in_bvhTree.nodes[currentNode.leftOffset + i], currentDepth + 1));
          }
        }
        else{
          if(currentDepth > bvhMaxDepth){
            bvhMaxDepth = currentDepth;
          }
          if(currentDepth < bvhMinDepth){
            bvhMinDepth = currentDepth;
          }
        }
      }

      if(pMinDepth != nullptr)
        (*pMinDepth) = bvhMinDepth;
      
      if(pMaxDepth != nullptr)
        (*pMaxDepth) = bvhMaxDepth;

      //std::cout << "\n[TestBVHdepth]:\nMininmum BVH depth:\t" << bvhMinDepth << "\nMaximum BVH depth:\t" << bvhMaxDepth <<std::endl;
      return (bvhMinDepth >= 1) && (bvhMaxDepth < 100);
    }
    
    /*
     * BVH Intervals test. Checks if the parent node actually points to the same indices.
     */
    static bool TestBVHintervals(const cbvh::BVHTree& in_bvhTree) {
      // If there's only one interval and therefore no child nodes, the test doesn't make much sense.
      if (in_bvhTree.intervals.size() == 1){
        //std::cout << "[TestBVHintervals]:\tPASSED! There's only one interval." << std::endl;
        return true;
      }

      for(int i = 0; i < in_bvhTree.intervals.size(); ++i)
      {
        int intervalStart  = in_bvhTree.intervals[i].start;
        int intervalLength = in_bvhTree.intervals[i].count;

        int accumulationLength = 0;
        
        if (cbvh::IsValid(in_bvhTree.nodes[i])) // if not leaf and not empty
        {
          for (int child_num = 0; child_num < 3; ++child_num) 
          {
            const auto childOffset = in_bvhTree.nodes[i].leftOffset + child_num;
            const auto nextIndex   = in_bvhTree.intervals[childOffset].start + in_bvhTree.intervals[childOffset].count;
            if (!cbvh::IsEmpty(in_bvhTree.nodes[childOffset + 1]) && nextIndex != in_bvhTree.intervals[childOffset + 1].start) 
            {
              //std::cout << "[TestBVHintervals]:\tFAILED! BVH intervals are calculated incorrectly." << std::endl;
              return false;
            }

            accumulationLength += in_bvhTree.intervals[in_bvhTree.nodes[i].leftOffset + child_num].count;
          }

          const auto accLen3 = (accumulationLength + in_bvhTree.intervals[in_bvhTree.nodes[i].leftOffset + 3].count);
          if (intervalLength != accLen3) 
          {
            //std::cout << "[TestBVHintervals](2):\t" << i << "\t" << intervalLength << " != " << accLen3 << std::endl;
            return false;
          }
        }
      
      }

      return true;
    }

    /*
     * BVH depth intervals test. Checks if the number of all nodes equals the number of nodes from 'depthRanges' depth levels.
     */
    static bool TestDepthIntervals(const cbvh::BVHTree& in_bvhTree) {
      int numberOfNodesByDepth = 0;
      for(int i = 0; i < in_bvhTree.depthRanges.size(); ++i){
        numberOfNodesByDepth += in_bvhTree.depthRanges[i].count;
      }
      auto nodeSize = in_bvhTree.nodes.size();
      return (nodeSize == numberOfNodesByDepth);
    }
    

    struct TraversalData
    {
      void FindLeafRec4(uint a_currOffser)
      {
        const cbvh::BVHNode& currNode = m_nodes[a_currOffser];
        
        if(IsLeaf(currNode))
        {
          const cbvh::Interval& currInterval = m_intervals[a_currOffser];
          for(int i=currInterval.start; i < currInterval.start + currInterval.count; i++)
          {
            m_triIndices.insert(m_indices[i*3+0]);
            m_triIndices.insert(m_indices[i*3+1]);
            m_triIndices.insert(m_indices[i*3+2]);
          }
        }
        else if(IsEmpty(currNode))
          return;
        else
        {
          FindLeafRec4(currNode.leftOffset + 0);
          FindLeafRec4(currNode.leftOffset + 1);
          FindLeafRec4(currNode.leftOffset + 2);
          FindLeafRec4(currNode.leftOffset + 3);
        }
      }

      void FindLeafRec2(uint a_currOffser)
      {
        const cbvh::BVHNode& currNode = m_nodes[a_currOffser];
        
        if(IsLeaf(currNode))
        {
          const cbvh::Interval& currInterval = m_intervals[a_currOffser];
          for(int i=currInterval.start; i < currInterval.start + currInterval.count; i++)
          {
            m_triIndices.insert(m_indices[i*3+0]);
            m_triIndices.insert(m_indices[i*3+1]);
            m_triIndices.insert(m_indices[i*3+2]);
          }
        }
        else if(IsEmpty(currNode))
          return;
        else
        {
          FindLeafRec2(currNode.leftOffset + 0);
          FindLeafRec2(currNode.leftOffset + 1);
        }
      }

      void FindLeafRec2Dyn(uint a_currOffser)
      {
        const cbvh::BVHNode& currNode = m_nodes[a_currOffser];
        
        if(IsLeaf(currNode))
        {
          const cbvh::Interval& currInterval = m_intervals[a_currOffser];
          for(int i=currInterval.start; i < currInterval.start + currInterval.count; i++)
          {
            m_triIndices.insert(m_indices[i*3+0]);
            m_triIndices.insert(m_indices[i*3+1]);
            m_triIndices.insert(m_indices[i*3+2]);
          }
        }
        else if(IsEmpty(currNode))
          return;
        else
        {
          FindLeafRec2Dyn(currNode.leftOffset  + 0);
          FindLeafRec2Dyn(currNode.escapeIndex + 0);
        }
      }

      int CountEmptyLeafes(uint a_currOffser)
      {
        const cbvh::BVHNode& currNode = m_nodes[a_currOffser];

        if(IsLeaf(currNode))
          return 0;
        else
        if(IsEmpty(currNode))
          return 1;
        else
          return CountEmptyLeafes(currNode.leftOffset + 0) + CountEmptyLeafes(currNode.leftOffset + 1) + 
                 CountEmptyLeafes(currNode.leftOffset + 2) + CountEmptyLeafes(currNode.leftOffset + 3);
      }

      bool EscapeIndexCheck(uint a_currOffset, uint a_parentEscapeIndex)
      {
        const cbvh::BVHNode& currNode = m_nodes[a_currOffset];
        
        if(!IsLeaf(currNode) && !IsEmpty(currNode))
        {
          const cbvh::BVHNode& child0 = m_nodes[currNode.leftOffset+0];
          const cbvh::BVHNode& child1 = m_nodes[currNode.leftOffset+1];
          const cbvh::BVHNode& child2 = m_nodes[currNode.leftOffset+2];
          const cbvh::BVHNode& child3 = m_nodes[currNode.leftOffset+3];

          const bool ckeck1 = (child0.escapeIndex == currNode.leftOffset+1);
          const bool ckeck2 = (child1.escapeIndex == currNode.leftOffset+2);
          const bool ckeck3 = (child2.escapeIndex == currNode.leftOffset+3);
          const bool ckeck4 = (child3.escapeIndex == a_parentEscapeIndex);

          const bool ckeck5 = EscapeIndexCheck(currNode.leftOffset+0, child0.escapeIndex);
          const bool ckeck6 = EscapeIndexCheck(currNode.leftOffset+1, child1.escapeIndex);
          const bool ckeck7 = EscapeIndexCheck(currNode.leftOffset+2, child2.escapeIndex);
          const bool ckeck8 = EscapeIndexCheck(currNode.leftOffset+3, child3.escapeIndex);

          return ckeck1 && ckeck2 && ckeck3 && ckeck4 && 
                 ckeck5 && ckeck6 && ckeck7 && ckeck8;
        }

        return true;
      }

      std::set<uint>        m_triIndices;
      const cbvh::BVHNode*  m_nodes;
      const cbvh::Interval* m_intervals;
      const uint*           m_indices;
    };

    static inline bool IsDegenerativeTriangle(size_t index, const uint* a_inputIndices, const float4* a_vertices)
    {
      Box4f box;
      const auto A = a_vertices[ a_inputIndices[index+0] ];
      const auto B = a_vertices[ a_inputIndices[index+1] ];
      const auto C = a_vertices[ a_inputIndices[index+2] ];
      box.include(A);
      box.include(B);
      box.include(C);
      const float4 boxSize = box.boxMax - box.boxMin; // #TODO: skip thin triangles also? make preset for this ??
      return (dot3f(boxSize,boxSize) <= 1e-24f) ;     // skip degenerative triangles
    }

    static inline bool IsLineTriangle(size_t index, const uint* a_inputIndices, const float4* a_vertices)
    {
      const auto A = a_vertices[ a_inputIndices[index+0] ];
      const auto B = a_vertices[ a_inputIndices[index+1] ];
      const auto C = a_vertices[ a_inputIndices[index+2] ];
      const float4 flatNorm = cross(A-B,A-C);
      return (dot3f(flatNorm,flatNorm) == 0.0f);
    }

    static bool TestLeafTrisSet(const cbvh::BVHTree& in_bvhTree, const uint* a_inputIndices, size_t a_inIndicesSize, const float4* a_vertices)
    {
      TraversalData data;
      data.m_nodes     = in_bvhTree.nodes.data();
      data.m_intervals = in_bvhTree.intervals.data();
      data.m_indices   = in_bvhTree.indicesReordered.data();
      data.m_triIndices.clear();
      
      if(in_bvhTree.format == cbvh::FMT_BVH4Node32_Interval32_Static)
        data.FindLeafRec4(0);
      else if(in_bvhTree.format == cbvh::FMT_BVH2Node32_Interval32_Static)
        data.FindLeafRec2(0);
      else if(in_bvhTree.format == cbvh::FMT_BVH2Node32_Interval32_Dynamic)
        data.FindLeafRec2Dyn(0);

      uint found1 = 0xFFFFFFFF;
      uint found2 = 0xFFFFFFFF;

      for(size_t i=0;i<a_inIndicesSize;i++)
      {
        uint index = a_inputIndices[i];
        
        if(data.m_triIndices.find(index) == data.m_triIndices.end())
        {
          bool point = IsDegenerativeTriangle(i,a_inputIndices,a_vertices); 
          bool line  = IsLineTriangle(i,a_inputIndices,a_vertices); 
          if(i%3==0 && !point && !line)
          {
             found1 = index;
             break;
          }
          else if(point)
          {
            i+=2;
            continue;
          }
        }
      }

      for(size_t i=0;i<in_bvhTree.indicesReordered.size();i++)
      {
        uint index = in_bvhTree.indicesReordered[i];
        
        if(data.m_triIndices.find(index) == data.m_triIndices.end())
        {
           found2 = index;
           break;
        }
      }

      // std::cout << "[TestLeafTrisSet]:\tFAILED!\t(" << found1 << ", " << found1 << ") " << std::endl;
      return (found1 == 0xFFFFFFFF && found2 == 0xFFFFFFFF);
    }

  static int TestEmptyNodes(const cbvh::BVHTree& in_bvhTree)
  {
    TraversalData data;
    data.m_nodes     = in_bvhTree.nodes.data();
    data.m_intervals = in_bvhTree.intervals.data();
    data.m_indices   = in_bvhTree.indicesReordered.data();
    data.m_triIndices.clear();
    return data.CountEmptyLeafes(0);
  }

  static bool TestEscapeIndex(const cbvh::BVHTree& in_bvhTree)
  {
    TraversalData data;
    data.m_nodes     = in_bvhTree.nodes.data();
    data.m_intervals = in_bvhTree.intervals.data();
    data.m_indices   = in_bvhTree.indicesReordered.data();
    data.m_triIndices.clear();
    return data.EscapeIndexCheck(0,-2);
  }

  void GenerateUniformRandomRays(const LiteMath::Box4f& a_scnBox, int a_numRays,
                                 LiteMath::float4* a_raysPos, LiteMath::float4* a_raysDir);


  bool TestRayTraversal(const cbvh::BVHTree& in_bvhTree, 
                        const LiteMath::float4* in_vPos4f, size_t a_maxVertexCount,
                        const uint* a_indicesOriginal, size_t a_indexNum);

  struct TraversalMetricsData
  {
    TraversalMetricsData() : avgNodesCount(0), avgLeafesCount(0), varLeafesCount(0), avgPrimsCount(0) {} 
    float avgNodesCount;
    float avgLeafesCount;
    float varLeafesCount; // variance, for LCV
    float avgPrimsCount;
  };

  TraversalMetricsData ComputeRayTraversalMetrics(const cbvh::BVHTree& in_bvhTree, 
                                                  const LiteMath::float4* in_vPos4f, size_t a_maxVertexCount,
                                                  const uint* a_indicesOriginal, size_t a_indexNum);

  struct TestSimpleMesh
  {
    static constexpr uint64_t POINTS_IN_TRIANGLE = 3;
    TestSimpleMesh(){}
    TestSimpleMesh(int a_vertNum, int a_indNum) { Resize(a_vertNum, a_indNum); }
    inline size_t VerticesNum()  const { return vPos4f.size(); }
    inline size_t IndicesNum()   const { return indices.size();  }
    inline size_t TrianglesNum() const { return IndicesNum() / POINTS_IN_TRIANGLE;  }
    inline void   Resize(uint32_t a_vertNum, uint32_t a_indNum) 
    {
      vPos4f.resize(a_vertNum);
      vNorm4f.resize(a_vertNum);
      vTang4f.resize(a_vertNum);
      vTexCoord2f.resize(a_vertNum);
      indices.resize(a_indNum);
      matIndices.resize(a_indNum/3); 
      assert(a_indNum%3 == 0); // PLEASE NOTE THAT CURRENT IMPLEMENTATION ASSUME ONLY TRIANGLE MESHES! 
    };
    
    std::vector<LiteMath::float4> vPos4f;      // 
    std::vector<LiteMath::float4> vNorm4f;     // 
    std::vector<LiteMath::float4> vTang4f;     // 
    std::vector<LiteMath::float2>  vTexCoord2f; // 
    std::vector<unsigned int>      indices;     // size = 3*TrianglesNum() for triangle mesh, 4*TrianglesNum() for quad mesh
    std::vector<unsigned int>      matIndices;  // size = 1*TrianglesNum()
  };

  std::vector<Box4f> TriangleBoxesFromMesh(const float4* a_vpos, size_t a_vertNum, const uint* a_indices, size_t a_indexNum);

  static void CheckBvhCorrectness_FMT2D(const cbvh::BVHTree& treeData, LiteMath::uint parent_index, LiteMath::uint child_index, int& stack_steps, const LiteMath::float4* vertices)
  {
    assert(treeData.format == cbvh::FMT_BVH2Node32_Interval32_Dynamic);

    LiteMath::float3 parentMin;
    LiteMath::float3 parentMax;

    assert(child_index != 0xffffffff);

    LiteMath::float3 childMin = treeData.nodes[child_index].boxMin;
    LiteMath::float3 childMax = treeData.nodes[child_index].boxMax;

    if (parent_index != 0xffffffff) 
    {
      assert(parent_index < treeData.nodes.size());
      parentMin = treeData.nodes[parent_index].boxMin;
      parentMax = treeData.nodes[parent_index].boxMax;

      assert((childMin.x >= parentMin.x) && (childMin.y >= parentMin.y) && (childMin.z >= parentMin.z));
      assert((childMax.x <= parentMax.x) && (childMax.y <= parentMax.y) && (childMax.z <= parentMax.z));
      stack_steps++;
    }

    if (treeData.nodes[child_index].leftOffset != 0xffffffff)  // check children
    {
      CheckBvhCorrectness_FMT2D(treeData, child_index, treeData.nodes[child_index].leftOffset, stack_steps, vertices);
      CheckBvhCorrectness_FMT2D(treeData, child_index, treeData.nodes[child_index].escapeIndex, stack_steps, vertices);
    } 
    else // check triangles into aabb 
    {
      uint32_t int_start = treeData.intervals[child_index].start;
      uint32_t int_end   = int_start + treeData.intervals[child_index].count;

      if (treeData.intervals[child_index].count < 1024)
      {
        for (uint32_t i = int_start; i < int_end; i++) 
        {
          uint32_t ind0 = treeData.indicesReordered[i * 3 + 0];
          uint32_t ind1 = treeData.indicesReordered[i * 3 + 1];
          uint32_t ind2 = treeData.indicesReordered[i * 3 + 2];

          const LiteMath::float4& tri_v0 = vertices[ind0];
          const LiteMath::float4& tri_v1 = vertices[ind1];
          const LiteMath::float4& tri_v2 = vertices[ind2];
            
          const bool inside0X = (tri_v0.x >= childMin.x) && (tri_v0.x <= childMax.x);
          const bool inside0Y = (tri_v0.y >= childMin.y) && (tri_v0.y <= childMax.y);
          const bool inside0Z = (tri_v0.z >= childMin.z) && (tri_v0.z <= childMax.z);
          
          const bool inside1X = (tri_v1.x >= childMin.x) && (tri_v1.x <= childMax.x);
          const bool inside1Y = (tri_v1.y >= childMin.y) && (tri_v1.y <= childMax.y);
          const bool inside1Z = (tri_v1.z >= childMin.z) && (tri_v1.z <= childMax.z);

          const bool inside2X = (tri_v2.x >= childMin.x) && (tri_v2.x <= childMax.x);
          const bool inside2Y = (tri_v2.y >= childMin.y) && (tri_v2.y <= childMax.y);
          const bool inside2Z = (tri_v2.z >= childMin.z) && (tri_v2.z <= childMax.z);

          assert(inside0X && inside0Y && inside0Z);
          assert(inside1X && inside1Y && inside1Z);
          assert(inside2X && inside2Y && inside2Z);
        }
      }
      else
      {
        std::cout << "ALERT! too many primitives in leaf: " <<  treeData.intervals[child_index].count << std::endl;
      }
    }

    return;
  }  

  static void CheckBvhCorrectness_FMT2S(const cbvh::BVHTree& treeData, LiteMath::uint parent_index, LiteMath::uint child_index, int& stack_steps, const LiteMath::float4* vertices)
  {
    assert(treeData.format == cbvh::FMT_BVH2Node32_Interval32_Static);

    LiteMath::float3 parentMin;
    LiteMath::float3 parentMax;

    assert(child_index != 0xffffffff);

    LiteMath::float3 childMin = treeData.nodes[child_index].boxMin;
    LiteMath::float3 childMax = treeData.nodes[child_index].boxMax;

    if (parent_index != 0xffffffff) 
    {
      assert(parent_index < treeData.nodes.size());
      parentMin = treeData.nodes[parent_index].boxMin;
      parentMax = treeData.nodes[parent_index].boxMax;

      assert((childMin.x >= parentMin.x) && (childMin.y >= parentMin.y) && (childMin.z >= parentMin.z));
      assert((childMax.x <= parentMax.x) && (childMax.y <= parentMax.y) && (childMax.z <= parentMax.z));
      stack_steps++;
    }

    if (treeData.nodes[child_index].leftOffset != 0xffffffff)  // check children
    {
      CheckBvhCorrectness_FMT2S(treeData, child_index, treeData.nodes[child_index].leftOffset+0, stack_steps, vertices);
      CheckBvhCorrectness_FMT2S(treeData, child_index, treeData.nodes[child_index].leftOffset+1, stack_steps, vertices);
    } 
    else // check triangles into aabb 
    {
      uint32_t int_start = treeData.intervals[child_index].start;
      uint32_t int_end   = int_start + treeData.intervals[child_index].count;

      if (treeData.intervals[child_index].count < 1024)
      {
        for (uint32_t i = int_start; i < int_end; i++) 
        {
          uint32_t ind0 = treeData.indicesReordered[i * 3 + 0];
          uint32_t ind1 = treeData.indicesReordered[i * 3 + 1];
          uint32_t ind2 = treeData.indicesReordered[i * 3 + 2];

          const LiteMath::float4& tri_v0 = vertices[ind0];
          const LiteMath::float4& tri_v1 = vertices[ind1];
          const LiteMath::float4& tri_v2 = vertices[ind2];
            
          const bool inside0X = (tri_v0.x >= childMin.x) && (tri_v0.x <= childMax.x);
          const bool inside0Y = (tri_v0.y >= childMin.y) && (tri_v0.y <= childMax.y);
          const bool inside0Z = (tri_v0.z >= childMin.z) && (tri_v0.z <= childMax.z);
          
          const bool inside1X = (tri_v1.x >= childMin.x) && (tri_v1.x <= childMax.x);
          const bool inside1Y = (tri_v1.y >= childMin.y) && (tri_v1.y <= childMax.y);
          const bool inside1Z = (tri_v1.z >= childMin.z) && (tri_v1.z <= childMax.z);

          const bool inside2X = (tri_v2.x >= childMin.x) && (tri_v2.x <= childMax.x);
          const bool inside2Y = (tri_v2.y >= childMin.y) && (tri_v2.y <= childMax.y);
          const bool inside2Z = (tri_v2.z >= childMin.z) && (tri_v2.z <= childMax.z);

          assert(inside0X && inside0Y && inside0Z);
          assert(inside1X && inside1Y && inside1Z);
          assert(inside2X && inside2Y && inside2Z);
        }
      }
      else
      {
        std::cout << "ALERT! too many primitives in leaf: " <<  treeData.intervals[child_index].count << std::endl;
      }
    }

    return;
  }  


  static void CheckBvhCorrectness_FMT4S(const cbvh::BVHTree& treeData, LiteMath::uint parent_index, LiteMath::uint child_index, int& stack_steps, const LiteMath::float4* vertices)
  {
    assert(treeData.format == cbvh::FMT_BVH4Node32_Interval32_Static);

    LiteMath::float3 parentMin;
    LiteMath::float3 parentMax;

    assert(child_index != 0xffffffff);

    LiteMath::float3 childMin = treeData.nodes[child_index].boxMin;
    LiteMath::float3 childMax = treeData.nodes[child_index].boxMax;

    if (parent_index != 0xffffffff) 
    {
      assert(parent_index < treeData.nodes.size());
      parentMin = treeData.nodes[parent_index].boxMin;
      parentMax = treeData.nodes[parent_index].boxMax;

      assert((childMin.x >= parentMin.x) && (childMin.y >= parentMin.y) && (childMin.z >= parentMin.z));
      assert((childMax.x <= parentMax.x) && (childMax.y <= parentMax.y) && (childMax.z <= parentMax.z));
      stack_steps++;
    }

    if (treeData.nodes[child_index].leftOffset < uint32_t(-3))  // check children
    {
      CheckBvhCorrectness_FMT4S(treeData, child_index, treeData.nodes[child_index].leftOffset+0, stack_steps, vertices);
      CheckBvhCorrectness_FMT4S(treeData, child_index, treeData.nodes[child_index].leftOffset+1, stack_steps, vertices);
      CheckBvhCorrectness_FMT4S(treeData, child_index, treeData.nodes[child_index].leftOffset+2, stack_steps, vertices);
      CheckBvhCorrectness_FMT4S(treeData, child_index, treeData.nodes[child_index].leftOffset+3, stack_steps, vertices);
    } 
    else // check triangles into aabb 
    {
      uint32_t int_start = treeData.intervals[child_index].start;
      uint32_t int_end   = int_start + treeData.intervals[child_index].count;

      if (treeData.intervals[child_index].count < 1024)
      {
        for (uint32_t i = int_start; i < int_end; i++) 
        {
          uint32_t ind0 = treeData.indicesReordered[i * 3 + 0];
          uint32_t ind1 = treeData.indicesReordered[i * 3 + 1];
          uint32_t ind2 = treeData.indicesReordered[i * 3 + 2];

          const LiteMath::float4& tri_v0 = vertices[ind0];
          const LiteMath::float4& tri_v1 = vertices[ind1];
          const LiteMath::float4& tri_v2 = vertices[ind2];
            
          const bool inside0X = (tri_v0.x >= childMin.x) && (tri_v0.x <= childMax.x);
          const bool inside0Y = (tri_v0.y >= childMin.y) && (tri_v0.y <= childMax.y);
          const bool inside0Z = (tri_v0.z >= childMin.z) && (tri_v0.z <= childMax.z);
          
          const bool inside1X = (tri_v1.x >= childMin.x) && (tri_v1.x <= childMax.x);
          const bool inside1Y = (tri_v1.y >= childMin.y) && (tri_v1.y <= childMax.y);
          const bool inside1Z = (tri_v1.z >= childMin.z) && (tri_v1.z <= childMax.z);

          const bool inside2X = (tri_v2.x >= childMin.x) && (tri_v2.x <= childMax.x);
          const bool inside2Y = (tri_v2.y >= childMin.y) && (tri_v2.y <= childMax.y);
          const bool inside2Z = (tri_v2.z >= childMin.z) && (tri_v2.z <= childMax.z);

          assert(inside0X && inside0Y && inside0Z);
          assert(inside1X && inside1Y && inside1Z);
          assert(inside2X && inside2Y && inside2Z);
        }
      }
      else
      {
        std::cout << "ALERT! too many primitives in leaf: " <<  treeData.intervals[child_index].count << std::endl;
      }
    }

    return;
  }  

  static bool TestTrisIndideBoxes(const cbvh::BVHTree& treeData, const float4* meshVertF4)
  {
    int stack_steps = 0;

    if(treeData.format == cbvh::FMT_BVH2Node32_Interval32_Dynamic)
      CheckBvhCorrectness_FMT2D(treeData, LiteMath::uint(0xffffffff), LiteMath::uint(0), stack_steps, meshVertF4);
    else if(treeData.format == cbvh::FMT_BVH2Node32_Interval32_Static)
      CheckBvhCorrectness_FMT2S(treeData, LiteMath::uint(0xffffffff), LiteMath::uint(0), stack_steps, meshVertF4);
    else if(treeData.format == cbvh::FMT_BVH4Node32_Interval32_Static)
      CheckBvhCorrectness_FMT4S(treeData, LiteMath::uint(0xffffffff), LiteMath::uint(0), stack_steps, meshVertF4);
    else
    {
      assert(false); // undefined tree format
    }

    return true;
  }

  static bool RunAllTests(cbvh::BVHTree& treeData, const LiteMath::float4* a_vpos, size_t a_vertNum, const uint* a_indices, size_t a_indexNum)
  {
    int depthMin = 0, depthMax = 0;
    bool test1 = cbvh_test::TestBVHdepth(treeData, &depthMin, &depthMax);
    if(test1)
      std::cout << "[TestBVHdepth]:\t\tPASSED!"; 
    else
      std::cout << "[TestBVHdepth]:\t\tFAILED!";
    std::cout << " (" << depthMin << "," << depthMax << ")" << std::endl;
    
    bool test2 = cbvh_test::TestBVHintervals(treeData);
    if(test2)
      std::cout << "[TestBVHintervals]:\tPASSED!" << std::endl;
    else
      std::cout << "[TestBVHintervals]:\tFAILED!" << std::endl;
    
    bool test3 = cbvh_test::TestDepthIntervals(treeData);
    if(test3)
      std::cout << "[TestDepthIntervals]:\tPASSED!" << std::endl;
    else
      std::cout << "[TestDepthIntervals]:\tFAILED!" << std::endl;
  
    bool test4 = cbvh_test::TestLeafTrisSet(treeData, a_indices, a_indexNum, a_vpos);
    if(test4)
      std::cout << "[TestLeafTrisSet]:\tPASSED!" << std::endl;
    else
      std::cout << "[TestLeafTrisSet]:\tFAILED!" << std::endl;
  
    int emptyNodes = cbvh_test::TestEmptyNodes(treeData);
    std::cout << "[TestEmptyNodes]:\tPASSED! / EmptyNodes = " << emptyNodes << std::endl;
  
    bool test5 = cbvh_test::TestEscapeIndex(treeData);
    if(test5)
      std::cout << "[TestEscapeIndex]:\tPASSED!" << std::endl;
    else
      std::cout << "[TestEscapeIndex]:\tFAILED!" << std::endl;
  
    bool test6 = cbvh_test::TestRayTraversal(treeData, a_vpos, a_vertNum, a_indices, a_indexNum);
    if(test6)
      std::cout << "[TestRayTraversal]:\tPASSED!" << std::endl;
    else
      std::cout << "[TestRayTraversal]:\tFAILED!" << std::endl;

    bool test7 = cbvh_test::TestTrisIndideBoxes(treeData, a_vpos);
    if(test7)
      std::cout << "[TestTrisIndideBoxes]:\tPASSED!" << std::endl;
    else
      std::cout << "[TestTrisIndideBoxes]:\tFAILED!" << std::endl;

    return test1 && test2 && test3 && test4 && test5 && test6 && test7;
  }
}