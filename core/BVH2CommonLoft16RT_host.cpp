#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>

#include "BVH2CommonLoft16RT.h"

#include "builders/cbvh_core.h"
#include "aligned_alloc.h"
using cbvh::BVHNode;
using cbvh::Interval;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr size_t reserveSize = 1000;

void BVH2_LOFT_16::ClearGeom()
{
  m_vertPos.reserve(std::max<size_t>(100000, m_vertPos.capacity()));
  m_indices.reserve(std::max<size_t>(100000 * 3, m_indices.capacity()));
  m_primIndices.reserve(std::max<size_t>(100000, m_primIndices.capacity()));

  m_vertPos.resize(0);
  m_indices.resize(0);
  m_primIndices.resize(0);

  m_allNodes.reserve(std::max<size_t>(100000, m_allNodes.capacity()));
  m_allNodes.resize(0);

  m_geomOffsets.reserve(std::max(reserveSize, m_geomOffsets.capacity()));
  m_geomOffsets.resize(0);

  m_geomBoxes.reserve(std::max<size_t>(reserveSize, m_geomBoxes.capacity()));
  m_geomBoxes.resize(0);

  m_bvhOffsets.reserve(std::max<size_t>(reserveSize, m_bvhOffsets.capacity()));
  m_bvhOffsets.resize(0);

  m_geomSize.resize(0);
  
  totalTrisMem = 0;
  totalTrisMem = 0;
  totalinstNumInst = 0;
  ClearScene();
}

size_t BVH2_LOFT_16::AppendTreeData(const std::vector<cbvh::BVHNode>& a_nodes, const std::vector<uint32_t>& a_indices, 
                                    const uint32_t *a_triIndices, size_t a_indNumber)
{
  const size_t oldSize = m_allNodes.size();
  size_t oldIndexSize  = m_indices.size();

  m_allNodes.insert(m_allNodes.end(), a_nodes.begin(), a_nodes.end());
  m_primIndices.insert(m_primIndices.end(), a_indices.begin(), a_indices.end());

  m_indices.resize(oldIndexSize + a_indices.size()*3);
  for(size_t i=0;i<a_indices.size();i++)
  {
    const uint32_t triId = a_indices[i];
    m_indices[oldIndexSize + 3*i+0] = a_triIndices[triId*3+0];
    m_indices[oldIndexSize + 3*i+1] = a_triIndices[triId*3+1];
    m_indices[oldIndexSize + 3*i+2] = a_triIndices[triId*3+2];
  }
  
  return oldSize;
}

uint32_t BVH2_LOFT_16::AddGeom_Triangles3f(const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride)
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

  // Build BVH for each geom and append it to big buffer
  //
  auto presets = cbvh2::BuilderPresetsFromString(m_builderName.c_str());
  auto bvhData = cbvh2::BuildBVH((const float*)(m_vertPos.data() + oldSizeVert), a_vertNumber, 16, a_triIndices, a_indNumber, presets);

  const size_t oldBvhSize = AppendTreeData(bvhData.nodes, bvhData.indices, a_triIndices, a_indNumber);
  m_bvhOffsets.push_back(uint32_t(oldBvhSize));

  totalTrisMem += (a_indNumber/3);
  m_geomSize.push_back(uint32_t(a_indNumber/3));
  return currGeomId;
}

void BVH2_LOFT_16::UpdateGeom_Triangles3f(uint32_t a_geomId, const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride)
{
  std::cout << "[BVH2_LOFT_16::UpdateGeom_Triangles3f]: "
            << "not implemeted!" << std::endl;
}

void BVH2_LOFT_16::ClearScene()
{
  m_instBoxes.reserve(std::max(reserveSize, m_instBoxes.capacity()));
  m_instMatricesInv.reserve(std::max(reserveSize, m_instMatricesInv.capacity()));
  m_instMatricesFwd.reserve(std::max(reserveSize, m_instMatricesFwd.capacity()));
  m_geomIdByInstId.reserve(std::max(reserveSize, m_geomIdByInstId.capacity()));

  m_instBoxes.resize(0);
  m_instMatricesInv.resize(0); // 
  m_instMatricesFwd.resize(0);
  m_geomIdByInstId.resize(0);

  if(m_tlasOffset != 0) // discard TLAS from 'm_allNodes' buffer
  { 
    m_allNodes.resize(m_tlasOffset);
    m_tlasOffset = 0;
  }
}

static inline half_float::half roundToNeg(float a_val) { return half_float::nextafter(half_float::half(a_val), -std::numeric_limits<half_float::half>::max()); }
static inline half_float::half roundToPos(float a_val) { return half_float::nextafter(half_float::half(a_val), +std::numeric_limits<half_float::half>::max()); }

static inline BVHNode16 ConvertNode(const BVHNode& a_node)
{
  BVHNode16 node;
  node.minXYZmaxX = half4(roundToNeg(a_node.boxMin.x), roundToNeg(a_node.boxMin.y), roundToNeg(a_node.boxMin.z), roundToPos(a_node.boxMax.x));
  node.maxYZ      = half2(roundToPos(a_node.boxMax.y), roundToPos(a_node.boxMax.z));
  node.leftOffset = a_node.leftOffset;
  return node; 
}

static inline std::vector<BVHNode16, aligned<BVHNode16, 64> > ConvertNodes(const std::vector<BVHNode, aligned<BVHNode, 64> >& a_nodes)
{
  std::vector<BVHNode16, aligned<BVHNode16, 64> > nodes16(a_nodes.size());
  for(size_t i=0;i<nodes16.size();i++)
    nodes16[i] = ConvertNode(a_nodes[i]);
  return nodes16;
}

void BVH2_LOFT_16::CommitScene(BuildQuality a_qualityLevel)
{
  m_tlasOffset = uint32_t(m_allNodes.size());

  cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_FAST, 1};
  //cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_MEDIUM, 1};
  auto nodesTLAS = cbvh2::BuildBVH((const cbvh::BVHNode*)m_instBoxes.data(), m_instBoxes.size(), presets);

  m_allNodes.insert(m_allNodes.end(), nodesTLAS.begin(), nodesTLAS.end());
  m_allNodes16 = ConvertNodes(m_allNodes);

  // reset stats
  //  
  m_stats.clear();
  m_stats.bvhTotalSize  = m_allNodes.size()*sizeof(BVHNode);
  m_stats.geomTotalSize = m_vertPos.size()*sizeof(float4) + m_indices.size()*sizeof(uint32_t);

  std::cout << "totalTrisMem  = " << totalTrisMem << std::endl;
  std::cout << "totalTrisInst = " << totalTrisInst << std::endl;
  std::cout << "totalNumInst  = " << totalinstNumInst << std::endl;
}

uint32_t BVH2_LOFT_16::AddInstance(uint32_t a_geomId, const float4x4 &a_matrix)
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
  
  totalinstNumInst++;
  totalTrisInst += m_geomSize[a_geomId];
  return oldSize;
}

void BVH2_LOFT_16::UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix)
{
  std::cout << "[BVH2_LOFT_16::UpdateInstance]: "
            << "not implemeted!" << std::endl;
}

ISceneObject *MakeBVH2CommonLoft16RT(const char *a_implName, const char* a_buildName) { return new BVH2_LOFT_16(a_buildName); }