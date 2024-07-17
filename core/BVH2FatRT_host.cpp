#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <memory>

#include "BVH2FatRT.h"
#include "nanort/nanort.h"
#include "builders/refitter.h"

using LiteMath::BBox3f;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr size_t reserveSize = 1000;

void BVH2FatRT::ClearGeom()
{
  m_vertPos.reserve(std::max<size_t>(100000, m_vertPos.capacity()));
  m_indices.reserve(std::max<size_t>(100000 * 3, m_indices.capacity()));
  m_primIndices.reserve(std::max<size_t>(100000, m_primIndices.capacity()));

  m_vertPos.resize(0);
  m_indices.resize(0);
  m_primIndices.resize(0);

  m_allNodesFat.reserve(std::max<size_t>(100000, m_allNodesFat.capacity()));
  m_allNodesFat.resize(0);

  m_geomOffsets.reserve(std::max(reserveSize, m_geomOffsets.capacity()));
  m_geomOffsets.resize(0);

  m_geomBoxes.reserve(std::max<size_t>(reserveSize, m_geomBoxes.capacity()));
  m_geomBoxes.resize(0);

  m_bvhOffsets.reserve(std::max<size_t>(reserveSize, m_bvhOffsets.capacity()));
  m_bvhOffsets.resize(0);

  ClearScene();
}

void BVH2FatRT::AppendTreeData(const std::vector<BVHNodeFat>& a_nodes, const std::vector<uint32_t>& a_indices, const uint32_t *a_triIndices, size_t a_indNumber)
{
  m_allNodesFat.insert(m_allNodesFat.end(), a_nodes.begin(), a_nodes.end());
  m_primIndices.insert(m_primIndices.end(), a_indices.begin(), a_indices.end());
  
  const size_t oldIndexSize  = m_indices.size();
  m_indices.resize(oldIndexSize + a_indices.size()*3);
  for(size_t i=0;i<a_indices.size();i++)
  {
    const uint32_t triId = a_indices[i];
    m_indices[oldIndexSize + 3*i+0] = a_triIndices[triId*3+0];
    m_indices[oldIndexSize + 3*i+1] = a_triIndices[triId*3+1];
    m_indices[oldIndexSize + 3*i+2] = a_triIndices[triId*3+2];
  }
}

uint32_t BVH2FatRT::AddGeom_Triangles3f(const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride)
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
  const size_t oldBvhSize = m_allNodesFat.size();
  m_bvhOffsets.push_back(oldBvhSize);
  
  auto presets = cbvh2::BuilderPresetsFromString(m_buildName.c_str());
  auto layout  = cbvh2::LayoutPresetsFromString(m_layoutName.c_str());
  auto bvhData = cbvh2::BuildBVHFat((const float*)(m_vertPos.data() + oldSizeVert), a_vertNumber, 16, a_triIndices, a_indNumber, presets, layout);
  

  //auto pRefitter = std::make_unique<Refitter>();

  AppendTreeData(bvhData.nodes, bvhData.indices, a_triIndices, a_indNumber);

  //pRefitter->Refit_BVH2Fat32(m_allNodesFat.data() + oldBvhSize, uint32_t(bvhData.nodes.size()), 
  //                           m_indices.data()     + oldSizeInd, uint32_t(bvhData.indices.size()), 
  //                           bvhData.depthRanges.data(), uint32_t(bvhData.depthRanges.size()),
  //                           (const float*)(m_vertPos.data() + oldSizeVert), a_vertNumber, 
  //                                          m_indices.data() + oldSizeInd, a_indNumber, vByteStride/4); // note that we must remember depth ranges for each object we are going to refit

  return currGeomId;
}

void BVH2FatRT::UpdateGeom_Triangles3f(uint32_t a_geomId, const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride)
{
  std::cout << "[BVH2FatRT::UpdateGeom_Triangles3f]: " << "not implemeted!" << std::endl; // not planned for this implementation (possible in general)
}

void BVH2FatRT::ClearScene()
{
  m_instBoxes.reserve(std::max(reserveSize, m_instBoxes.capacity()));
  m_instMatricesInv.reserve(std::max(reserveSize, m_instMatricesInv.capacity()));
  m_instMatricesFwd.reserve(std::max(reserveSize, m_instMatricesFwd.capacity()));

  m_geomIdByInstId.reserve(std::max(reserveSize, m_geomIdByInstId.capacity()));

  m_instBoxes.resize(0);
  m_instMatricesInv.resize(0);
  m_instMatricesFwd.resize(0);
  m_geomIdByInstId.resize(0);

  m_firstSceneCommit = true;
}

void DebugPrintNodes(const std::vector<BVHNode>& nodes, const std::string& a_fileName)
{
  std::ofstream fout(a_fileName.c_str());

  for(size_t i=0;i<nodes.size();i++)
  {
    const auto& currBox = nodes[i];
    fout << "node[" << i << "]:" << std::endl;
    fout << "  bmin = { " << currBox.boxMin[0] << " " << currBox.boxMin[1] << " " << currBox.boxMin[2] << " } | " << currBox.leftOffset  << std::endl;
    fout << "  bmax = { " << currBox.boxMax[0] << " " << currBox.boxMax[1] << " " << currBox.boxMax[2] << " } | " << currBox.escapeIndex << std::endl;
  } 
}

void DebugPrintBoxes(const std::vector<Box4f>& nodes, const std::string& a_fileName)
{
  std::ofstream fout(a_fileName.c_str());

  for(size_t i=0;i<nodes.size();i++)
  {
    const auto& currBox = nodes[i];
    fout << "node[" << i << "]:" << std::endl;
    fout << "  bmin = { " << currBox.boxMin[0] << " " << currBox.boxMin[1] << " " << currBox.boxMin[2] << " " << currBox.boxMin[3]  << std::endl;
    fout << "  bmax = { " << currBox.boxMax[0] << " " << currBox.boxMax[1] << " " << currBox.boxMax[2] << " " << currBox.boxMax[3] << std::endl;
  } 
}

void BVH2FatRT::CommitScene(BuildQuality a_qualityLevel)
{
  if(m_firstSceneCommit)
  //if(true)
  {
    //cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_RIGHT, cbvh2::BVH_CONSTRUCT_FAST, 1};
    //cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_ROPES, cbvh2::BVH_CONSTRUCT_FAST, 1};
    //cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_MEDIUM, 1};
    cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_FAST, 1};
    m_nodesTLAS = cbvh2::BuildBVH((const cbvh2::BVHNode*)m_instBoxes.data(), m_instBoxes.size(), presets);
    
    //DebugPrintNodes(m_nodesTLAS, "z01_tlas.txt");
    //DebugPrintBoxes(m_instBoxes, "y01_boxes.txt");

    // reset stats
    //  
    m_stats.clear();
    m_stats.bvhTotalSize  = m_allNodesFat.size()*sizeof(BVHNodeFat) + m_nodesTLAS.size()*sizeof(BVHNode);
    m_stats.geomTotalSize = m_vertPos.size()*sizeof(float4)         + m_indices.size()*sizeof(uint32_t);
    m_firstSceneCommit    = false;

    m_pLBVHBuilder->Reserve(m_instBoxes.size(), cbvh2::BVH2_LEFT_OFFSET);
    m_pLBVHBuilder->CommitDeviceData(); // not needed for CPU in fact
    m_tlasReordered.resize(m_instBoxes.size());
  }
  else
  {
    // update 'm_nodesTLAS' from 'm_instBoxes' without additional allocations; 
    // currentlu use 'm_tlasReordered' as dummy buffer (according to single BLAS per TLAS node restriction)
    //
    m_pLBVHBuilder->BuildFromBoxes((const float4*)(m_instBoxes.data()), uint32_t(m_instBoxes.size()), 
                                    m_nodesTLAS.data(), uint32_t(m_nodesTLAS.size()), m_tlasReordered.data());

    //DebugPrintNodes(m_nodesTLAS, "z02_tlas.txt");
    //DebugPrintBoxes(m_instBoxes, "y02_boxes.txt");                                   
  }
}

uint32_t BVH2FatRT::AddInstance(uint32_t a_geomId, const float4x4 &a_matrix)
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

void BVH2FatRT::UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix)
{
  if(a_instanceId > m_geomIdByInstId.size())
  {
    std::cout << "[BVH2FatRT::UpdateInstance]: " << "bad instance id == " << a_instanceId << "; size == " << m_geomIdByInstId.size() << std::endl;
    return;
  }

  const uint32_t geomId = m_geomIdByInstId[a_instanceId];
  const float4 boxMin   = m_geomBoxes[geomId].boxMin;
  const float4 boxMax   = m_geomBoxes[geomId].boxMax;

  // (1) mult mesh bounding box vertices with matrix to form new bouding box for instance
  float4 boxVertices[8]{
      a_matrix * float4{boxMin.x, boxMin.y, boxMin.z, 1.0f},
      a_matrix * float4{boxMax.x, boxMin.y, boxMin.z, 1.0f},
      a_matrix * float4{boxMin.x, boxMax.y, boxMin.z, 1.0f},
      a_matrix * float4{boxMin.x, boxMin.y, boxMax.z, 1.0f},
      a_matrix * float4{boxMax.x, boxMax.y, boxMin.z, 1.0f},
      a_matrix * float4{boxMax.x, boxMin.y, boxMax.z, 1.0f},
      a_matrix * float4{boxMin.x, boxMax.y, boxMax.z, 1.0f},
      a_matrix * float4{boxMax.x, boxMax.y, boxMax.z, 1.0f},
  };

  Box4f newBox;
  for (size_t i = 0; i < 8; i++)
    newBox.include(boxVertices[i]);

  m_instBoxes      [a_instanceId] = newBox;
  m_instMatricesFwd[a_instanceId] = a_matrix;
  m_instMatricesInv[a_instanceId] = inverse4x4(a_matrix);
}

ISceneObject* MakeBVH2FatRT(const char* a_implName, const char* a_buildName, const char* a_layoutName) 
{
  return new BVH2FatRT(a_buildName, a_layoutName); 
}
