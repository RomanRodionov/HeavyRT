#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <chrono>
#include <cctype>

#include "BVH2FatCompactRT.h"

constexpr size_t reserveSize = 1000;

void BVH2FatRTCompact::ClearGeom()
{
  m_vertPos.reserve(std::max<size_t>(100000, m_vertPos.capacity()));
  m_indices.reserve(std::max<size_t>(100000 * 3, m_indices.capacity()));
  m_primIndices.reserve(std::max<size_t>(100000, m_primIndices.capacity()));

  m_vertPos.resize(0);
  m_indices.resize(0);
  m_primIndices.resize(0);

  m_allNodesFat16B.reserve(std::max<size_t>(100000, m_allNodesFat16B.capacity()));
  m_allNodesFat16B.resize(0);

  m_geomOffsets.reserve(std::max(reserveSize, m_geomOffsets.capacity()));
  m_geomOffsets.resize(0);

  m_geomBoxes.reserve(std::max<size_t>(reserveSize, m_geomBoxes.capacity()));
  m_geomBoxes.resize(0);

  m_bvhOffsets.reserve(std::max<size_t>(reserveSize, m_bvhOffsets.capacity()));
  m_bvhOffsets.resize(0);

  ClearScene();
}

uint32_t BVH2FatRTCompact::AddGeom_Triangles3f(const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride)
{
  const size_t vStride = vByteStride / 4;
  assert(vByteStride % 4 == 0);

  const uint32_t currGeomId = uint32_t(m_geomOffsets.size());
  const size_t oldSizeVert  = m_vertPos.size();
  const size_t oldSizeInd   = m_indices.size();

  m_geomOffsets.push_back(uint2(oldSizeInd, oldSizeVert));
  m_vertPos.resize(oldSizeVert + a_vertNumber);

  Box4f bbox;
  for (size_t i = 0; i < a_vertNumber; i++)
  {
    const float4 v = float4(a_vpos3f[i * vStride + 0], a_vpos3f[i * vStride + 1], a_vpos3f[i * vStride + 2], 1.0f);
    m_vertPos[oldSizeVert + i] = v;
    bbox.include(v);
  }

  m_geomBoxes.push_back(bbox);

  // Build BVH for each geom and append it to big buffer;
  // append data to global arrays and fix offsets
  //
  if(m_layoutName != "SuperTreeletAlignedMerged2" && m_layoutName != "SuperTreeletAlignedMerged4" && m_layoutName != "SuperTreeletAlignedMerged8")
  {
    std::cout << "[BVH2FatRTCompact]::WARNING: Bad layout for CALBVH: " << m_layoutName.c_str() << std::endl;
    m_layoutName = "SuperTreeletAlignedMerged4";
  }

  auto presets = cbvh2::BuilderPresetsFromString(m_buildName.c_str());
  auto layout  = cbvh2::LayoutPresetsFromString(m_layoutName.c_str());
  auto bvhData = cbvh2::BuildBVHFatCompressed((const float*)(m_vertPos.data() + oldSizeVert), a_vertNumber, 16, a_triIndices, a_indNumber, presets, layout);
  
  const uint oldBvhSize  = m_allNodesFat16B.size();
  m_bvhOffsets.push_back(oldBvhSize);

  if(bvhData.nodes.size() % (layout.grSzXX*2) != 0)
    std::cout << "[BVH2FatRTCompact]: ALERT! Bad size of CALBVH node array: bvhData.nodes.size() = " << bvhData.nodes.size() << std::endl;

  // append tree and geomatry data
  {
    m_log2_group_size = bvhData.log2_group_size;
    m_allNodesFat16B.insert(m_allNodesFat16B.end(), bvhData.nodes.begin(), bvhData.nodes.end());
    m_primIndices.insert(m_primIndices.end(), bvhData.indices.begin(), bvhData.indices.end());
    
    size_t oldIndexSize  = m_indices.size();
    m_indices.resize(oldIndexSize + bvhData.indices.size()*3);
    for(size_t i=0;i<bvhData.indices.size();i++)
    {
      const uint32_t triId = bvhData.indices[i];
      m_indices[oldIndexSize + 3*i+0] = a_triIndices[triId*3+0];
      m_indices[oldIndexSize + 3*i+1] = a_triIndices[triId*3+1];
      m_indices[oldIndexSize + 3*i+2] = a_triIndices[triId*3+2];
    }
  }

  return currGeomId;
}

void BVH2FatRTCompact::UpdateGeom_Triangles3f(uint32_t a_geomId, const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride)
{
  std::cout << "[BVH2FatRTCompact::UpdateGeom_Triangles3f]: "
            << "not implemeted!" << std::endl;
}

uint32_t BVH2FatRTCompact::AddInstance(uint32_t a_geomId, const float4x4 &a_matrix)
{
  const auto &box = m_geomBoxes[a_geomId];

  // (1) mult mesh bounding box vertices with matrix to form new bouding box for instance
  float4 boxVertices[8]{
      a_matrix * float4{box.boxMin.x, box.boxMin.y, box.boxMin.z, 1.0f},

      a_matrix * float4{box.boxMax.x, box.boxMin.y, box.boxMin.z, 1.0f},
      a_matrix * float4{box.boxMin.x, box.boxMax.y, box.boxMin.z, 1.0f},
      a_matrix * float4{box.boxMin.x, box.boxMin.y, box.boxMax.z, 1.0f},

      a_matrix * float4{box.boxMax.x, box.boxMax.y, box.boxMin.z, 1.0f},
      a_matrix * float4{box.boxMax.x, box.boxMin.y, box.boxMax.z, 1.0f},
      a_matrix * float4{box.boxMin.x, box.boxMax.y, box.boxMax.z, 1.0f},

      a_matrix * float4{box.boxMax.x, box.boxMax.y, box.boxMax.z, 1.0f},
  };

  Box4f newBox;
  for (size_t i = 0; i < 8; i++)
    newBox.include(boxVertices[i]);

  // (2) append bounding box and matrices
  //
  const uint32_t oldSize = uint32_t(m_instBoxes.size());

  m_instBoxes.push_back(newBox);
  m_instMatricesFwd.push_back(a_matrix);
  m_instMatricesInv.push_back(inverse4x4(a_matrix));
  m_geomIdByInstId.push_back(a_geomId);

  return oldSize;
}

void BVH2FatRTCompact::UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix)
{
  std::cout << "[BVH2FatRTCompact::UpdateInstance]: "
            << "not implemeted!" << std::endl;
}

void BVH2FatRTCompact::ClearScene()
{
  m_instBoxes.reserve(std::max(reserveSize, m_instBoxes.capacity()));
  m_instMatricesInv.reserve(std::max(reserveSize, m_instMatricesInv.capacity()));
  m_instMatricesFwd.reserve(std::max(reserveSize, m_instMatricesFwd.capacity()));
  m_geomIdByInstId.reserve(std::max(reserveSize, m_geomIdByInstId.capacity()));

  m_instBoxes.resize(0);
  m_instMatricesInv.resize(0);
  m_instMatricesFwd.resize(0);
  m_geomIdByInstId.resize(0);
}

void BVH2FatRTCompact::CommitScene(BuildQuality a_qualityLevel)
{ 
  //cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_RIGHT, cbvh2::BVH_CONSTRUCT_FAST, 1};
  //cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_ROPES, cbvh2::BVH_CONSTRUCT_FAST, 1};
  cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_MEDIUM, 1};
  m_nodesTLAS = cbvh2::BuildBVH((const cbvh2::BVHNode*)m_instBoxes.data(), m_instBoxes.size(), presets);
  
  m_stats.clear();
  m_stats.bvhTotalSize  = m_allNodesFat16B.size()*sizeof(uint4) + m_nodesTLAS.size()*sizeof(BVHNode);
  m_stats.geomTotalSize = m_vertPos.size()*sizeof(float4) + m_indices.size()*sizeof(uint32_t);
}

ISceneObject* MakeBVH2FatRTCompact(const char* a_implName, const char* a_buildName, const char* a_layoutName)
{
  return new BVH2FatRTCompact(a_buildName, a_layoutName);
}
