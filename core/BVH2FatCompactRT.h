#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <chrono>

#define HALFFLOAT
#include "BVH2FatRT.h"

using cbvh2::BVHNode;
using cbvh2::BVHNodeFat;

struct BVH2FatRTCompact : public ISceneObject //: public BVH2FatRT
{
  BVH2FatRTCompact(const char* a_buildName, const char* a_layoutName) : m_buildName (a_buildName  != nullptr ? a_buildName : ""), 
                                                                        m_layoutName(a_layoutName != nullptr ? a_layoutName : "") {}
  ~BVH2FatRTCompact() override {}

  const char* Name() const override { return "BVH2FatCompact"; }
  const char* BuildName() const override { return m_buildName.c_str(); };

  void ClearGeom() override;

  void CommitScene(BuildQuality a_qualityLevel) override;

  uint32_t AddGeom_Triangles3f(const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride) override;
  void     UpdateGeom_Triangles3f(uint32_t a_geomId, const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride) override;

  void     ClearScene() override;
  uint32_t AddInstance(uint32_t a_geomId, const float4x4 &a_matrix) override;
  void     UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix) override;

  CRT_Hit RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar) override;
  bool    RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar) override;

  uint32_t GetGeomNum() const override { return uint32_t(m_geomBoxes.size()); }
  uint32_t GetInstNum() const override { return uint32_t(m_instBoxes.size()); }
  const LiteMath::float4* GetGeomBoxes() const override { return (const LiteMath::float4*)m_geomBoxes.data(); }

//protected:

  void BVH2TraverseF32(const float3 ray_pos, const float3 ray_dir, float tNear,
                       uint32_t instId, uint32_t geomId, uint32_t stack[STACK_SIZE], bool stopOnFirstHit,
                       CRT_Hit* pHit);

  void IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
                                    float tNear, uint32_t instId, uint32_t geomId,
                                    uint32_t a_start, uint32_t a_count,
                                    CRT_Hit *pHit);

  std::vector<uint4> m_allNodesFat16B;
  uint               m_log2_group_size;
  
  std::vector<Box4f> m_geomBoxes;
  std::vector<Box4f> m_instBoxes;

  std::vector<float4x4> m_instMatricesInv; ///< inverse instance matrices
  std::vector<float4x4> m_instMatricesFwd; ///< instance matrices

  std::vector<float4>     m_vertPos;
  std::vector<uint32_t>   m_indices;
  std::vector<uint32_t>   m_primIndices;

  std::vector<BVHNode>    m_nodesTLAS;
  std::vector<uint32_t>   m_bvhOffsets;

  std::vector<uint2>      m_geomOffsets;
  std::vector<uint32_t>   m_geomIdByInstId;

  // Format name the tree build from
  std::string m_buildName;
  std::string m_layoutName;
};

ISceneObject* MakeBVH2FatRTCompact(const char* a_implName, const char* a_buildName, const char* a_layoutName);

