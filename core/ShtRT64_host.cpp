#include "ShtRT64.h"
#include <numeric>
#include <map>

constexpr size_t reserveSize = 1000;

void ShtRT64::ClearGeom()
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

  ClearScene();
}

void ShtRT64::AppendTreeData(const std::vector<BVHNode> &a_nodes, const std::vector<uint32_t> &a_indices,
                           const uint32_t *a_triIndices, size_t a_indNumber)
{
  m_allNodes.insert(m_allNodes.end(), a_nodes.begin(), a_nodes.end());
  m_primIndices.insert(m_primIndices.end(), a_indices.begin(), a_indices.end());

  std::vector<uint32_t> indicesReordered(a_indices.size() * 3); // #TODO: opt this, may directly append to 'm_indices'
  for(uint32_t i = 0; i < uint32_t(a_indices.size()); i++)
  {
    const uint32_t triId = a_indices[i];
    indicesReordered[3 * i + 0] = a_triIndices[triId * 3 + 0];
    indicesReordered[3 * i + 1] = a_triIndices[triId * 3 + 1];
    indicesReordered[3 * i + 2] = a_triIndices[triId * 3 + 2];
  }
  m_indices.insert(m_indices.end(), indicesReordered.begin(), indicesReordered.end());
}

uint32_t
ShtRT64::AddGeom_Triangles3f(const float *a_vpos3f, size_t a_vertNumber, const uint32_t *a_triIndices, size_t a_indNumber,
                           BuildQuality a_qualityLevel, size_t vByteStride)
{
  const size_t vStride = vByteStride / 4;
  assert(vByteStride % 4 == 0);

  const uint32_t currGeomId = uint32_t(m_geomOffsets.size());
  const size_t oldSizeVert = m_vertPos.size();
  const size_t oldSizeInd = m_indices.size();

  m_geomOffsets.push_back(uint2(oldSizeInd, oldSizeVert));

  m_vertPos.resize(oldSizeVert + a_vertNumber);

  Box4f bbox;
  for(size_t i = 0; i < a_vertNumber; i++)
  {
    const float4 v = float4(a_vpos3f[i * vStride + 0], a_vpos3f[i * vStride + 1], a_vpos3f[i * vStride + 2], 1.0f);
    m_vertPos[oldSizeVert + i] = v;
    bbox.include(v);
  }

  m_geomBoxes.push_back(bbox);

  // Build BVH for each geom and append it to big buffer
  //
  cbvh::BVHPresets presets = GetShtPresetsFromName(m_builderName, a_qualityLevel);

  // append data to global arrays and fix offsets
  //
  const size_t oldBvhSize = m_allNodes.size(); // both Nodes and Intervals always have same size, we know that
  m_bvhOffsets.push_back(oldBvhSize);

  //auto currentOffset = m_bvhOffsets[m_bvhOffsets.size() - 1];
  auto bvhData = cbvh::BuildBVH(m_vertPos.data() + oldSizeVert, a_vertNumber, a_triIndices, a_indNumber, presets);
  BuildBVHWithUncles(bvhData.nodes, 0, 0, 0);
  AppendTreeData(bvhData.nodes, bvhData.indicesReordered, a_triIndices, a_indNumber);

  // convert tree to C32 ('Compact32' format)
  //
  for(size_t i = 0; i < bvhData.intervals.size(); i++)
    if(int(m_allNodes[oldBvhSize + i].leftOffset) < 0)
      m_allNodes[oldBvhSize + i].leftOffset = PackOffsetAndSize(bvhData.intervals[i].start,
                                                                bvhData.intervals[i].count);


  return currGeomId;
}

void ShtRT64::UpdateGeom_Triangles3f(uint32_t a_geomId, const float *a_vpos3f, size_t a_vertNumber,
                                   const uint32_t *a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel,
                                   size_t vByteStride)
{
  std::cout << "[BVH2ShtRT64::UpdateGeom_Triangles3f]: "
            << "not implemeted!" << std::endl;
}

void ShtRT64::ClearScene()
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

void ShtRT64::CommitScene(BuildQuality a_qualityLevel)
{
  // Init tables


  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  InitHashTables();

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms for InitHashTables"
            << std::endl;

//    visualize_tree(m_allNodes);

  // Build LBVH
  //
  cbvh2::BuilderPresets presets = {cbvh2::BVH2_LEFT_OFFSET, cbvh2::BVH_CONSTRUCT_MEDIUM, 1};
  m_nodesTLAS = cbvh2::BuildBVH((const cbvh::BVHNode*)m_instBoxes.data(), m_instBoxes.size(), presets);

  //Count stats
  #ifndef KERNEL_SLICER
  m_stats.clear();
  m_stats.bvhTotalSize  = m_allNodes.size() * sizeof(BVHNode) + m_nodesTLAS.size()*sizeof(BVHNode) + HashMemSize();
  m_stats.geomTotalSize = m_vertPos.size()  * sizeof(float4)  + m_indices.size()*sizeof(uint32_t);
  #endif
}

uint32_t ShtRT64::AddInstance(uint32_t a_geomId, const float4x4 &a_matrix)
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
  for(size_t i = 0; i < 8; i++)
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


void ShtRT64::UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix)
{
  std::cout << "[ShtRT64::UpdateInstance]: "
            << "not implemeted!" << std::endl;
}

ISceneObject *MakeShtRT64(const char *a_implName, const char *a_buildName)
{ return new ShtRT64(a_buildName); }

cbvh::BVHPresets ShtRT64::GetShtPresetsFromName(const std::string &a_buildName, BuildQuality a_qualityLevel)
{
  cbvh::BVHPresets presets;
  presets.childrenNum = 2;
  presets.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Static;
  presets.btype = (a_qualityLevel == BuildQuality::BUILD_LOW) ? cbvh::BVH_CONSTRUCT_FAST : cbvh::BVH_CONSTRUCT_MEDIUM;
  presets.useEmbreeIfPossiable = false;
  if(a_buildName.find("cbvh_hq") != std::string::npos)
    presets.btype = cbvh::BVH_CONSTRUCT_QUALITY;
  if(a_buildName.find("cbvh_hqs") != std::string::npos || a_buildName.find("cbvh_meds") != std::string::npos)
    presets.enableSpatialSplit = true;
  if(a_buildName.find("cbvh_embree") != std::string::npos)
    presets.btype = cbvh::BVH_CONSTRUCT_EMBREE2;
  char symb = a_buildName[a_buildName.size() - 1];
  if(std::isdigit(symb))
    presets.primsInLeaf = int(symb) - int('0');
  return presets;
}

void
ShtRT64::BuildBVHWithUncles(std::vector<BVHNode> &bvhNodes, uint32_t node, uint32_t uncle, uint32_t brother)
{
  bvhNodes[node].escapeIndex = uncle;
  auto leftOffset  = bvhNodes[node].leftOffset;
  auto rightOffset = leftOffset + 1;

  if(leftOffset != EMPTY_NODE && !(leftOffset & LEAF_BIT))
  {
    BuildBVHWithUncles(bvhNodes, leftOffset, brother, rightOffset);
    if(rightOffset != EMPTY_NODE && !(rightOffset & LEAF_BIT))
    {
      BuildBVHWithUncles(bvhNodes, rightOffset, brother, leftOffset);
    }
  }
}

// Remove hash tables
ShtRT64::~ShtRT64()
{
#ifdef LINEAR_HASHING
  hashTable.resize(0);
  displacementTables.resize(0);
#else
  for(size_t i = 0; i < hashTables.size(); ++i)
  {
    delete hashTables[i];
  }
  #endif
}

#ifdef LINEAR_HASHING

void ShtRT64::InitHashTables()
{
  const size_t objNumber = m_bvhOffsets.size();

  std::cout << "Starting init hash tables..." << std::endl;

  std::vector<uint32_t> sizes;

  for (uint32_t i = 0; i < objNumber - 1; ++i) {
    sizes.push_back(m_bvhOffsets[i + 1] - m_bvhOffsets[i]);
  }
  dispSizes.resize(objNumber);
  hSizes.resize(objNumber);
  sizes.push_back(m_allNodes.size() - m_bvhOffsets[objNumber - 1]);
  dOffsets.push_back(0);
  hOffsets.push_back(0);
  for (uint32_t objectId = 0; objectId < objNumber; ++objectId) {
    uint32_t numKeys = 0;
    uint32_t offset = m_bvhOffsets[objectId];
    uint32_t size = sizes[objectId];
    std::map<uint64_t, uint32_t> keysAddressMap;

    // count existing node keys and associated addresses
    CountNodesKeys(offset, 0, 1, numKeys, keysAddressMap);

    // Displacement size is numOfKeys // 2 and lower to 2^n
    CountDisplacementSize(2*size, objectId);

    // Hash size is lower co_prime number for numOfKey * 2
    CountHashtableSize(2*size, objectId);
    auto dSize = dispSizes[objectId];
    auto hSize = hSizes[objectId];
    displacementTables.resize(displacementTables.size() + dSize);
    hashTable.resize(hashTable.size() + hSize);
    dOffsets.push_back(displacementTables.size());
    hOffsets.push_back(hashTable.size());
    std::map<uint32_t, std::vector<uint64_t>> displacementMap;

    // We count entries: {dispValue : listOf(keys)} such that key % D == dispValue
    InitDisplacementMap(displacementMap, keysAddressMap, dSize);

    std::vector<std::pair<uint32_t, std::vector<uint64_t>>> displacementVector;
    for (auto &dispPair: displacementMap) {
      displacementVector.push_back(dispPair);
    }

    displacementMap.erase(displacementMap.begin(), displacementMap.end());

    InitHashTable(displacementVector, keysAddressMap, objectId);
  }
}

void ShtRT64::InitDisplacementMap(std::map<uint32_t, std::vector<uint64_t>> &displacementMap, std::map<uint64_t, uint32_t> &keys, uint32_t dSize) {
  for (auto &keyPair: keys) {
    auto key = keyPair.first;
    auto value = key & (dSize - 1);
    if (displacementMap.find(value) != displacementMap.end()) {
      displacementMap[value].push_back(key);
    }
    else {
      displacementMap[value] = {key};
    }
  }
}

bool ShtRT64::CanPlaceValuesToHashtable(std::vector<uint64_t> &keys, uint32_t dispValue, uint32_t hSize, uint32_t hOffset) {
  for (auto &key: keys) {
    auto hash = (key + dispValue) % hSize;
    if (hashTable[hOffset + hash] != 0) {
      return false;
    }
  }
  return true;
}
#include <set>

void ShtRT64::InitHashTable(std::vector<std::pair<uint32_t, std::vector<uint64_t>>> &displacementVector, std::map<uint64_t, uint32_t> &keysAddressMap, uint32_t objectId) {
  uint32_t hSize = hSizes[objectId];
  uint32_t hOffset = hOffsets[objectId];
  uint32_t dOffset = dOffsets[objectId];
  std::set<uint32_t> hashes;
  for (auto &i: displacementVector) {
    uint32_t dispValue = 0;
    uint32_t dispIndex = i.first;
    auto &keys = i.second;
    std::set<uint32_t> localProcessed;
    while(!CanPlaceValuesToHashtable(keys, dispValue, hSize, hOffset)) {
      dispValue++;
    }
    for (auto &key: keys) {
      displacementTables[dOffset + dispIndex] = dispValue;
      auto hash = (key + dispValue) % hSize;
      hashes.insert(hash);
      hashTable[hOffset + hash] = keysAddressMap[key];
    }
  }
}

void ShtRT64::CountDisplacementSize(uint32_t counter, uint32_t objectId) {
  auto v = counter / 2;
  v -= 1;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v += 1;
  while (v >= counter / 2) {
      v /= 2;
  }
  if (v == 0) {
    v = 1;
  }
  dispSizes[objectId] = v;
}

void ShtRT64::CountHashtableSize(uint32_t counter, uint32_t objectId) {
  auto v = counter * 2;
  auto currentValue = v + 1;
  while(true) {
      if (std::gcd(currentValue, v) == 1) {
          hSizes[objectId] = currentValue;
          break;
      }
      currentValue += 1;
  }
}

void ShtRT64::CountNodesKeys(uint32_t offset, const uint32_t currentNode, const uint64_t nodeKey, uint32_t &counter, std::map<uint64_t, uint32_t> &keysAddressMap) {
  keysAddressMap[nodeKey] = currentNode;
  counter++;
  const auto leftNode = m_allNodes[offset + currentNode].leftOffset;
  const auto rightNode = m_allNodes[offset + currentNode].leftOffset + 1;
  if(leftNode != EMPTY_NODE && !(leftNode & LEAF_BIT))
  {
    CountNodesKeys(offset, EXTRACT_OFFSET(leftNode), nodeKey << 1, counter, keysAddressMap);
    if(rightNode != EMPTY_NODE && !(rightNode & LEAF_BIT))
    {
      CountNodesKeys(offset, EXTRACT_OFFSET(rightNode), nodeKey << 1 ^ 1, counter, keysAddressMap);
    }
  }
}
#else
void ShtRT64::InitHashTables() {
    const size_t nodesNumber = m_allNodes.size();
    const size_t objNumber = m_bvhOffsets.size();

    std::vector<size_t> sizes;
    sizes.resize(objNumber);
    for (size_t i = 0; i < (objNumber - 1); ++i) {
        sizes[i] = m_bvhOffsets[i + 1] - m_bvhOffsets[i];
    }
    sizes[objNumber - 1] = nodesNumber - m_bvhOffsets[objNumber - 1];

    hashTables.resize(objNumber);

    //init hash table
    for (size_t i = 0; i < objNumber; ++i) {

        hashTables[i] = new MockHashTable(sizes[i], m_bvhOffsets[i], m_allNodes);
    }
}
#endif
