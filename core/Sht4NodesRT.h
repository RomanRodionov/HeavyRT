#ifndef LOADER_APP_SHT4NODESRT_H
#define LOADER_APP_SHT4NODESRT_H
#include <iostream>
#include "CrossRT.h"
#include "raytrace_common.h"
#include "builders/cbvh_core.h"
#include <map>
#include <chrono>

#define LINEAR_HASHING

using LiteMath::float4x4;
using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;
using cbvh::BVHNode;
using LiteMath::Box4f;
using LiteMath::uint2;
using LiteMath::uint4;


class Sht4NodesRT : public ISceneObject
{
public:
  Sht4NodesRT(const char *a_builderName = "nanort") : m_builderName(a_builderName)
  {
    
  };

  ~Sht4NodesRT() override;

  [[nodiscard]] const char *Name() const override
  { return "Sht4NodesRT"; }

  void ClearGeom() override;

  uint32_t
  AddGeom_Triangles3f(const float *a_vpos3f,
                      size_t a_vertNumber,
                      const uint32_t *a_triIndices,
                      size_t a_indNumber,
                      BuildQuality a_qualityLevel,
                      size_t vByteStride) override;

  void
  UpdateGeom_Triangles3f(uint32_t a_geomId,
                         const float *a_vpos3f,
                         size_t a_vertNumber,
                         const uint32_t *a_triIndices,
                         size_t a_indNumber,
                         BuildQuality a_qualityLevel,
                         size_t vByteStride) override;

  void ClearScene() override;

  void CommitScene(BuildQuality a_qualityLevel) override;

  uint32_t AddInstance(uint32_t a_geomId,
                       const float4x4 &a_matrix) override;

  void UpdateInstance(uint32_t a_instanceId,
                      const float4x4 &a_matrix) override;

  CRT_Hit RayQuery_NearestHit(float4 posAndNear,
                              float4 dirAndFar) override;

  bool RayQuery_AnyHit(float4 posAndNear,
                       float4 dirAndFar) override;

  uint32_t GetGeomNum() const override { return uint32_t(m_geomBoxes.size()); }
  uint32_t GetInstNum() const override { return uint32_t(m_instBoxes.size()); }
  const LiteMath::float4* GetGeomBoxes() const override { return (const LiteMath::float4*)m_geomBoxes.data(); }

//protected:

  std::vector<Box4f> m_geomBoxes;
  std::vector<Box4f> m_instBoxes;

  std::vector<float4x4> m_instMatricesInv; ///< inverse instance matrices
  std::vector<float4x4> m_instMatricesFwd; ///< instance matrices

  std::vector<float4> m_vertPos; /// vertices positions
  std::vector<uint32_t> m_indices; /// triangle indices
  std::vector<uint32_t> m_primIndices;


  std::vector<BVHNode> m_nodesTLAS;
  std::vector<BVHNode> m_allNodes;
  std::vector<uint32_t> m_bvhOffsets;

  std::vector<uint2> m_geomOffsets;
  std::vector<uint32_t> m_geomIdByInstId;
  
  size_t HashMemSize() const { return sizeof(uint32_t)*(hashTable.size() + displacementTables.size() + hSizes.size() + dispSizes.size() + hOffsets.size() + dOffsets.size()); }
  std::vector<uint32_t> hashTable;
  std::vector<uint32_t> displacementTables;
  std::vector<uint32_t> hSizes;
  std::vector<uint32_t> dispSizes;
  std::vector<uint32_t> hOffsets;
  std::vector<uint32_t> dOffsets;

  const char *m_builderName;

  static uint32_t CountTrailingZeros(uint64_t trail);

  void IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
                                    float tNear, uint32_t instId, uint32_t geomId,
                                    uint32_t a_start, uint32_t a_count,
                                    CRT_Hit *pHit);

  void Traverse(const float3 rayPos, const float3 ray_dir, float tMin,
                uint32_t instId, uint32_t geomId, CRT_Hit *hit);

  void InitHashTables();

  void InitHashTable(std::vector<std::pair<uint32_t, std::vector<uint64_t>>> &displacementVector, std::map<uint64_t, uint32_t> &keysAddressMap, uint32_t objectId);

  void InitDisplacementMap(std::map<uint32_t, std::vector<uint64_t>> &displacementMap, std::map<uint64_t, uint32_t> &keys, uint32_t dSize);

  void CountNodesKeys(uint32_t offset, const uint32_t currentNode, const uint64_t nodeKey, uint32_t &counter, std::map<uint64_t, uint32_t> &keysAddressMap);

  void CountDisplacementSize(uint32_t size, uint32_t objectId);

  bool CanPlaceValuesToHashtable(std::vector<uint64_t> &keys, uint32_t dispValue, uint32_t hSize, uint32_t hOffset);

  void CountHashtableSize (uint32_t size, uint32_t objectId);

  virtual void AppendTreeData(const std::vector<BVHNode> &a_nodes,
                              const std::vector<uint32_t> &a_indices,
                              const uint32_t *a_triIndices, size_t a_indNumber);
};


#endif //LOADER_APP_SHT4NODESRT_H
