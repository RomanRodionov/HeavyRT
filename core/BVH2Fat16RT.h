#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>

#define HALFFLOAT
#include "LiteMath.h"
#include "aligned_alloc.h"

using LiteMath::cross;
using LiteMath::dot;
using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;
using LiteMath::float4x4;
using LiteMath::int2;
using LiteMath::inverse4x4;
using LiteMath::normalize;
using LiteMath::sign;
using LiteMath::to_float3;
using LiteMath::uint;
using LiteMath::uint2;
using LiteMath::uint4;
using LiteMath::Box4f;

#include "CrossRT.h"
#include "raytrace_common.h"
#include "builders/cbvh.h"

using cbvh2::BVHNode;
using cbvh2::BVHNodeFat;
using cbvh2::BVHNodeFat16;

// main class
//
struct BVH2Fat16RT : public ISceneObject
{
  BVH2Fat16RT(const char* a_buildName, const char* a_layoutName) : m_buildName(a_buildName != nullptr ? a_buildName : ""), 
                                                                 m_layoutName(a_layoutName != nullptr ? a_layoutName : "") {}
  ~BVH2Fat16RT() override {}

  const char* Name() const override { return "BVH2Fat16RT"; }
  const char* BuildName() const override { return m_buildName.c_str(); };

  void ClearGeom() override;

  virtual uint32_t AddGeom_Triangles3f(const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride) override;
  void UpdateGeom_Triangles3f(uint32_t a_geomId, const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride) override;

  void ClearScene() override;
  virtual void CommitScene(BuildQuality a_qualityLevel) override;

  uint32_t AddInstance(uint32_t a_geomId, const float4x4 &a_matrix) override;
  void UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix) override;

  CRT_Hit RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar) override;
  bool RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar) override;

  uint32_t GetGeomNum() const override { return uint32_t(m_geomBoxes.size()); }
  uint32_t GetInstNum() const override { return uint32_t(m_instBoxes.size()); }
  const LiteMath::float4* GetGeomBoxes() const override { return (const LiteMath::float4*)m_geomBoxes.data(); }

//protected:

  void IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
                                    float tNear, uint32_t instId, uint32_t geomId,
                                    uint32_t a_start, uint32_t a_count,
                                    CRT_Hit *pHit);

  virtual void BVH2TraverseF16(const float3 ray_pos, const float3 ray_dir, float tNear, 
                               uint32_t instId, uint32_t geomId, uint32_t stack[STACK_SIZE], bool stopOnFirstHit,
                               CRT_Hit *pHit);

                     
  std::vector<Box4f> m_geomBoxes;
  std::vector<Box4f> m_instBoxes;

  std::vector<float4x4> m_instMatricesInv; ///< inverse instance matrices
  std::vector<float4x4> m_instMatricesFwd; ///< instance matrices

  std::vector<float4>   m_vertPos;
  std::vector<uint32_t> m_indices;
  std::vector<uint32_t> m_primIndices;

  std::vector<BVHNode>    m_nodesTLAS;
  std::vector<BVHNodeFat16, aligned<BVHNodeFat16, 32> > m_allNodesFat16;
  std::vector<uint32_t>   m_bvhOffsets;

  std::vector<uint2>    m_geomOffsets;
  std::vector<uint32_t> m_geomIdByInstId;

  // Format name the tree build from
  const std::string m_buildName;
  const std::string m_layoutName;

};
