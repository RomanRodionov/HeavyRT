#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>

#include "CrossRT.h"
#include "raytrace_common.h"

#include "builders/cbvh.h"

using LiteMath::dot;
using LiteMath::sign;
using LiteMath::cross;
using LiteMath::float4x4;
using LiteMath::uint2;
using LiteMath::uint4;
using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;
using LiteMath::normalize;
using LiteMath::inverse4x4;
using LiteMath::to_float3;
using LiteMath::as_float;
using LiteMath::as_uint;

using LiteMath::Box4f;
using cbvh2::BVHNode;

struct BVH4CommonRT : public ISceneObject
{
  BVH4CommonRT(const char* a_builderName = "cbvh_embree2") :  m_builderName(a_builderName) {}
  ~BVH4CommonRT() override {}
  
  const char* Name() const override { return "BVH4Common"; }

  void ClearGeom() override;

  uint32_t AddGeom_Triangles3f(const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride) override;
  void     UpdateGeom_Triangles3f(uint32_t a_geomId, const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride) override;

  void ClearScene() override; 
  void CommitScene(BuildQuality a_qualityLevel) override; 
  
  uint32_t AddInstance(uint32_t a_geomId, const float4x4& a_matrix) override;
  void     UpdateInstance(uint32_t a_instanceId, const float4x4& a_matrix) override;

  CRT_Hit  RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar) override;
  bool     RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar) override;
  
  uint32_t GetGeomNum() const override { return uint32_t(m_geomBoxes.size()); }
  uint32_t GetInstNum() const override { return uint32_t(m_instBoxes.size()); }
  const LiteMath::float4* GetGeomBoxes() const override { return (const LiteMath::float4*)m_geomBoxes.data(); }
  
//protected:
  
  void IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir, 
                                    float tNear, uint32_t instId, uint32_t geomId,
                                    uint32_t a_start, uint32_t a_count,
                                    CRT_Hit* pHit);

  void BVH4Traverse(const float3 ray_pos, const float3 ray_dir, float tNear, uint32_t instId, uint32_t geomId, uint32_t stack[STACK_SIZE],
                    bool stopOnFirstHit, CRT_Hit* pHit);

  void AppendTreeData(const std::vector<BVHNode>& a_nodes, const std::vector<uint32_t>& a_indices, 
                      const uint32_t *a_triIndices, size_t a_indNumber);

  std::vector<Box4f>    m_geomBoxes;
  std::vector<Box4f>    m_instBoxes;

  std::vector<float4x4> m_instMatricesInv; ///< inverse instance matrices
  std::vector<float4x4> m_instMatricesFwd; ///< instance matrices

  std::vector<float4>   m_vertPos;
  std::vector<uint32_t> m_indices;
  std::vector<uint32_t> m_primIndices;

  std::vector<BVHNode>  m_allNodes;
  std::vector<uint32_t> m_bvhOffsets;
  std::vector<BVHNode>  m_nodesTLAS;

  //std::vector<uint2>     m_vertStartSize;
  //std::vector<uint2>     m_indStartSize;
  std::vector<uint32_t>  m_geomIdByInstId;
  std::vector<uint2>     m_geomOffsets;

  std::string m_builderName;
};
