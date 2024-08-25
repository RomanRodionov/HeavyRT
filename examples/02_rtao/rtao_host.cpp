#include <cfloat>
#include "rtao.h"
#include "../loader/hydraxml.h"
#include "../loader/cmesh.h"
#include "loader/gltf_loader.h"

#include <chrono>
#include "qmc_sobol_niederreiter.h"

using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;

using LiteMath::perspectiveMatrix;
using LiteMath::lookAt;
using LiteMath::inverse4x4;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

RTAO::RTAO() {}

void RTAO::Render(uint32_t* a_outColor, uint32_t a_width, uint32_t a_height, const char* a_what, int a_passNum)
{
  CalcAOBlock(a_outColor, a_width, a_height, a_passNum);
}

void RTAO::SetAORadius(float radius) {
  m_aoMaxRadius = m_presets.isAORadiusInMeters ? radius : std::abs(m_aoBoxSize)*radius;
}

#if defined(__ANDROID__)
bool RTAO::LoadScene(const char* a_scenePath, AAssetManager* assetManager)
#else
bool RTAO::LoadScene(const char* a_scenePath)
#endif
{
  m_pAccelStruct->ClearGeom();

  const std::string& path = a_scenePath;

  auto isHydraScene = gltf_loader::ends_with(path, ".xml");

  if(isHydraScene)
  {
#if defined(__ANDROID__)
    return LoadSceneHydra(path, assetManager);
#else
    return LoadSceneHydra(path);
#endif
  }
  else
  {
#if defined(__ANDROID__)
    return LoadSceneGLTF(path, assetManager);
#else
    return LoadSceneGLTF(path);
#endif
  }

  return false;
}


#ifdef __ANDROID__
bool RTAO::LoadSceneHydra(const std::string& a_path, AAssetManager* assetManager)
#else
bool RTAO::LoadSceneHydra(const std::string& a_path)
#endif
{
  hydra_xml::HydraScene scene;
  if(
#if defined(__ANDROID__)
scene.LoadState(assetManager, a_path) < 0

#else
    scene.LoadState(a_path) < 0
#endif
          )
    return false;

  for(auto cam : scene.Cameras())
  {
    float aspect   = float(m_width) / float(m_height);
    auto proj      = perspectiveMatrix(cam.fov, aspect, cam.nearPlane, cam.farPlane);
    auto worldView = lookAt(float3(cam.pos), float3(cam.lookAt), float3(cam.up));

    m_zNearFar     = float2(cam.nearPlane, cam.farPlane);
    m_projInv      = inverse4x4(proj);
    m_worldViewInv = inverse4x4(worldView);

    break; // take first cam
  }

  m_vertOffset.reserve(1024);
  m_vNorm4f.resize(0);
  m_vPos4f.resize(0);

  std::vector<LiteMath::Box4f> meshBoxes;
  meshBoxes.reserve(100);
  
  m_totalTris = 0;
  std::vector<uint64_t> trisPerObject;
  trisPerObject.reserve(1000);

  m_pAccelStruct->ClearGeom();
  for(auto meshPath : scene.MeshFiles())
  {
    std::cout << "[LoadScene]: mesh = " << meshPath.c_str() << std::endl;
#if defined(__ANDROID__)
    auto currMesh = cmesh::LoadMeshFromVSGF2(assetManager, meshPath.c_str());
#else
    auto currMesh = cmesh::LoadMeshFromVSGF2(meshPath.c_str());
#endif
    auto geomId   = m_pAccelStruct->AddGeom_Triangles3f((const float*)currMesh.vPos4f.data(), currMesh.vPos4f.size(),
                                                        currMesh.indices.data(), currMesh.indices.size(), BUILD_HIGH, sizeof(float)*4);
    (void)geomId; // silence "unused variable" compiler warnings
    m_totalTris += currMesh.indices.size()/3;
    trisPerObject.push_back(currMesh.indices.size()/3);
    // we need this to estimate mesh bounding boxes
    {
      const float4* vPos4f = (const float4*)currMesh.vPos4f.data();
      LiteMath::Box4f box;
      for(size_t i=0;i<currMesh.vPos4f.size();i++)
        box.include(vPos4f[i]);
      meshBoxes.push_back(box);
    }

    m_matIdOffsets.push_back(m_matIdByPrimId.size());
    m_vertOffset.push_back(m_vNorm4f.size());

    m_matIdByPrimId.insert(m_matIdByPrimId.end(), currMesh.matIndices.begin(), currMesh.matIndices.end());
    m_triIndices.insert(m_triIndices.end(), currMesh.indices.begin(), currMesh.indices.end());

    m_vNorm4f.insert(m_vNorm4f.end(), currMesh.vNorm4f.begin(), currMesh.vNorm4f.end());
    m_vPos4f.insert(m_vPos4f.end(), currMesh.vPos4f.begin(), currMesh.vPos4f.end());
  }
  
  m_totalTrisVisiable = 0;
  Box4f sceneBox;
  m_normMatrices.clear(); m_normMatrices.reserve(1000);
  m_pAccelStruct->ClearScene();
  for(auto inst : scene.InstancesGeom())
  {
    m_pAccelStruct->AddInstance(inst.geomId, inst.matrix);
    m_normMatrices.push_back(transpose(inverse4x4(inst.matrix)));
    m_totalTrisVisiable += trisPerObject[inst.geomId];
    {
      const auto &box = meshBoxes[inst.geomId];
      // (1) mult mesh bounding box vertices with matrix to form new bouding box for instance
      float4 boxVertices[8]{
          inst.matrix * float4{box.boxMin.x, box.boxMin.y, box.boxMin.z, 1.0f},
          inst.matrix * float4{box.boxMax.x, box.boxMin.y, box.boxMin.z, 1.0f},
          inst.matrix * float4{box.boxMin.x, box.boxMax.y, box.boxMin.z, 1.0f},
          inst.matrix * float4{box.boxMin.x, box.boxMin.y, box.boxMax.z, 1.0f},
          inst.matrix * float4{box.boxMax.x, box.boxMax.y, box.boxMin.z, 1.0f},
          inst.matrix * float4{box.boxMax.x, box.boxMin.y, box.boxMax.z, 1.0f},
          inst.matrix * float4{box.boxMin.x, box.boxMax.y, box.boxMax.z, 1.0f},
          inst.matrix * float4{box.boxMax.x, box.boxMax.y, box.boxMax.z, 1.0f},
      };
      for (size_t i = 0; i < 8; i++)
         sceneBox.include(boxVertices[i]);
    }
 
  } 
  m_pAccelStruct->CommitScene();

  float4 boxSize = sceneBox.boxMax - sceneBox.boxMin;
  m_aoBoxSize    = std::max(boxSize.x, std::max(boxSize.y, boxSize.z));
  SetAORadius(m_presets.aoRayLength);
  m_aoRaysCount  = m_presets.aoRaysNum;

  InitQMCTable();

  std::cout << "[RTAO::LoadScene] finish qmc table init.   " << std::endl;
  return true;
}

struct DataRefs
{
  DataRefs(std::unordered_map<int, std::pair<uint32_t, BBox3f>> &a_loadedMeshesToMeshId, std::vector<uint64_t>& a_trisPerObject, Box4f& a_sceneBBox,
           int& a_gltfCam, float4x4& a_worldViewInv, std::shared_ptr<ISceneObject> a_pAccelStruct, uint64_t& a_totalTris, uint64_t& a_totalTrisVisible,
           std::vector<uint32_t>& a_matIdOffsets, std::vector<uint32_t>& a_vertOffset, std::vector<uint32_t>& a_matIdByPrimId,
           std::vector<uint32_t>& a_triIndices, std::vector<float4>&   a_vNorm4f, std::vector<float4>&   a_vPos4f, std::vector<LiteMath::float4x4>& a_normMatrices) :
          m_loadedMeshesToMeshId(a_loadedMeshesToMeshId), m_trisPerObject(a_trisPerObject), m_sceneBBox(a_sceneBBox),
          m_gltfCamId(a_gltfCam), m_worldViewInv(a_worldViewInv), m_pAccelStruct(a_pAccelStruct), m_totalTris(a_totalTris), m_totalTrisVisible(a_totalTrisVisible),
          m_matIdOffsets(a_matIdOffsets), m_vertOffset(a_vertOffset), m_matIdByPrimId(a_matIdByPrimId),
          m_triIndices(a_triIndices), m_vNorm4f(a_vNorm4f), m_vPos4f(a_vPos4f), m_normMatrices(a_normMatrices)
  {}
  std::unordered_map<int, std::pair<uint32_t, BBox3f>>& m_loadedMeshesToMeshId;
  std::vector<uint64_t>&                                m_trisPerObject;
  Box4f&                                                m_sceneBBox;

  int&                          m_gltfCamId;
  float4x4&                     m_worldViewInv;
  std::shared_ptr<ISceneObject> m_pAccelStruct;
  uint64_t&                     m_totalTris;
  uint64_t&                     m_totalTrisVisible;

  std::vector<uint32_t>& m_matIdOffsets;
  std::vector<uint32_t>& m_vertOffset;
  std::vector<uint32_t>& m_matIdByPrimId;
  std::vector<uint32_t>& m_triIndices;
  std::vector<float4>&   m_vNorm4f;
  std::vector<float4>&   m_vPos4f;
  std::vector<LiteMath::float4x4>& m_normMatrices;
};

static void LoadGLTFNodesRecursive(const tinygltf::Model &a_model, const tinygltf::Node& a_node, const LiteMath::float4x4& a_parentMatrix,
                                   const DataRefs& a_refs)
{
  auto nodeMatrix = a_parentMatrix * gltf_loader::transformMatrixFromGLTFNode(a_node);

  for (size_t i = 0; i < a_node.children.size(); i++)
  {
    LoadGLTFNodesRecursive(a_model, a_model.nodes[a_node.children[i]], nodeMatrix, a_refs);
  }

  if(a_node.camera > -1 && a_node.camera == a_refs.m_gltfCamId)
  {
    // works only for simple cases ?
    float3 eye = {0, 0, 0};
    float3 center = {0, 0, -1};
    float3 up = {0, 1, 0};
    auto tmp = lookAt(nodeMatrix * eye, nodeMatrix * center, up);
    a_refs.m_worldViewInv = inverse4x4(tmp);
  }
  if(a_node.mesh > -1)
  {
    if(!a_refs.m_loadedMeshesToMeshId.count(a_node.mesh))
    {
      const tinygltf::Mesh mesh = a_model.meshes[a_node.mesh];
      auto simpleMesh = gltf_loader::simpleMeshFromGLTFMesh(a_model, mesh);

      if(simpleMesh.VerticesNum() > 0)
      {
        auto meshId = a_refs.m_pAccelStruct->AddGeom_Triangles3f((const float*)simpleMesh.vPos4f.data(), simpleMesh.vPos4f.size(),
                                                                 simpleMesh.indices.data(), simpleMesh.indices.size(),
                                                                 BUILD_HIGH, sizeof(float)*4);
        a_refs.m_loadedMeshesToMeshId[a_node.mesh] = {meshId, simpleMesh.bbox};

        std::cout << "Loading mesh # " << meshId << std::endl;

        a_refs.m_totalTris += simpleMesh.indices.size() / 3;
        a_refs.m_trisPerObject.push_back(simpleMesh.indices.size() / 3);

        a_refs.m_matIdOffsets.push_back(a_refs.m_matIdByPrimId.size());
        a_refs.m_vertOffset.push_back(a_refs.m_vNorm4f.size());

        a_refs.m_matIdByPrimId.insert(a_refs.m_matIdByPrimId.end(), simpleMesh.matIndices.begin(), simpleMesh.matIndices.end() );
        a_refs.m_triIndices.insert(a_refs.m_triIndices.end(), simpleMesh.indices.begin(), simpleMesh.indices.end());

        a_refs.m_vNorm4f.insert(a_refs.m_vNorm4f.end(), simpleMesh.vNorm4f.begin(),     simpleMesh.vNorm4f.end());
        a_refs.m_vPos4f.insert(a_refs.m_vPos4f.end(), simpleMesh.vPos4f.begin(), simpleMesh.vPos4f.end());
      }
    }
    a_refs.m_totalTrisVisible += a_refs.m_trisPerObject[a_refs.m_loadedMeshesToMeshId[a_node.mesh].first];
    auto tmp_box = a_refs.m_loadedMeshesToMeshId[a_node.mesh].second;
    Box4f mesh_box {};
    mesh_box.boxMin = to_float4(tmp_box.boxMin, 1.0);
    mesh_box.boxMax = to_float4(tmp_box.boxMax, 1.0);
    Box4f inst_box = {};
    inst_box.boxMin = nodeMatrix * mesh_box.boxMin;
    inst_box.boxMax = nodeMatrix * mesh_box.boxMax;

    a_refs.m_sceneBBox.include(inst_box);
    a_refs.m_pAccelStruct->AddInstance(a_refs.m_loadedMeshesToMeshId[a_node.mesh].first, nodeMatrix);
    a_refs.m_normMatrices.push_back(transpose(inverse4x4(nodeMatrix)));
  }
}


#ifdef __ANDROID__
bool RTAO::LoadSceneGLTF(const std::string& a_path, AAssetManager* assetManager)
#else
bool RTAO::LoadSceneGLTF(const std::string& a_path)
#endif
{
  tinygltf::Model gltfModel;
  tinygltf::TinyGLTF gltfContext;
  std::string error, warning;

#ifdef __ANDROID__
  tinygltf::asset_manager = assetManager;
#endif

  bool loaded = gltfContext.LoadASCIIFromFile(&gltfModel, &error, &warning, a_path);

  if(!loaded)
  {
    std::cerr << "Cannot load glTF scene from: " << a_path << std::endl;;
    return false;
  }

  const tinygltf::Scene& scene = gltfModel.scenes[0];

  float aspect   = float(m_width) / float(m_height);
  int gltfCamId = -1;
  for(size_t i = 0; i < gltfModel.cameras.size(); ++i)
  {
    auto gltfCam = gltfModel.cameras[i];
    if(gltfCam.type == "perspective")
    {
      auto proj = perspectiveMatrix(gltfCam.perspective.yfov / DEG_TO_RAD, aspect,
                                    gltfCam.perspective.znear, gltfCam.perspective.zfar);
      m_zNearFar = float2(gltfCam.perspective.znear, gltfCam.perspective.zfar);
      m_projInv = inverse4x4(proj);

      gltfCamId = i;
      break; // take first perspective cam
    }
    else if (gltfCam.type == "orthographic")
    {
      std::cerr << "Orthographic camera not supported!" << std::endl;
    }
  }

  // load and instance geometry
  {
    std::vector<uint64_t> trisPerObject;
    trisPerObject.reserve(1000);
    m_totalTrisVisiable = 0;

    m_matIdOffsets.reserve(1024);
    m_vertOffset.reserve(1024);
    m_matIdByPrimId.reserve(128000);
    m_triIndices.reserve(128000*3);

    m_vNorm4f.resize(0);
    m_vPos4f.resize(0);

    m_totalTris = 0;
    m_pAccelStruct->ClearGeom();

    m_normMatrices.clear();
    m_normMatrices.reserve(1000);
    m_pAccelStruct->ClearScene();

    Box4f sceneBBox;
    std::unordered_map<int, std::pair<uint32_t, BBox3f>> loaded_meshes_to_meshId;
    for(size_t i = 0; i < scene.nodes.size(); ++i)
    {
      const tinygltf::Node node = gltfModel.nodes[scene.nodes[i]];
      auto identity = LiteMath::float4x4();
      //LoadGLTFNodesRecursive(gltfModel, node, identity, loaded_meshes_to_meshId, trisPerObject, sceneBBox);
      LoadGLTFNodesRecursive(gltfModel, node, identity,
                             DataRefs(loaded_meshes_to_meshId, trisPerObject, sceneBBox,
                                      gltfCamId, m_worldViewInv, m_pAccelStruct, m_totalTris, m_totalTrisVisiable,
                                      m_matIdOffsets, m_vertOffset, m_matIdByPrimId, m_triIndices,
                                      m_vNorm4f, m_vPos4f, m_normMatrices));
    }

    // glTF scene can have no cameras specified
    if(gltfCamId == -1)
    {
      std::tie(m_worldViewInv, m_projInv) = gltf_loader::makeCameraFromSceneBBox(m_width, m_height, sceneBBox);
    }
    
    const float3 camPos    = m_worldViewInv*float3(0,0,0);
    const float4 camLookAt = LiteMath::normalize(m_worldViewInv.get_col(2));
    std::cout << "[GLTF]: camPos2    = (" << camPos.x << "," << camPos.y << "," << camPos.z << ")" << std::endl;
    std::cout << "[GLTF]: camLookAt2 = (" << camLookAt.x << "," << camLookAt.y << "," << camLookAt.z << ")" << std::endl;
    std::cout << "[GLTF]: scnBoxMin  = (" << sceneBBox.boxMin.x << "," << sceneBBox.boxMin.y << "," << sceneBBox.boxMin.z << ")" << std::endl;
    std::cout << "[GLTF]: scnBoxMax  = (" << sceneBBox.boxMax.x << "," << sceneBBox.boxMax.y << "," << sceneBBox.boxMax.z << ")" << std::endl;

    m_pAccelStruct->CommitScene();

    float4 boxSize = sceneBBox.boxMax - sceneBBox.boxMin;
    m_aoBoxSize    = std::max(boxSize.x, std::max(boxSize.y, boxSize.z));
    SetAORadius(m_presets.aoRayLength);
    m_aoRaysCount  = m_presets.aoRaysNum;

    InitQMCTable();

    std::cout << "[RTAO::LoadScene] finish qmc table init.   " << std::endl;
  }

  return true;
}


void RTAO::InitQMCTable()
{
  std::cout << "[RTAO::LoadScene] start  qmc table init... " << std::endl;
  m_aoRandomsTile.resize(AO_TILE_SIZE*AO_TILE_SIZE*m_aoRaysCount);
  std::vector<uint32_t> counters(AO_TILE_SIZE*AO_TILE_SIZE*m_aoRaysCount);
  for(auto& cointer : counters)
    cointer = 0;
  
  unsigned int table[QRNG_DIMENSIONS][QRNG_RESOLUTION];
  initQuasirandomGenerator(table);

  for(size_t i=0;i<m_aoRandomsTile.size();i++)
  {
    const float x = rndQmcSobolN(uint32_t(i), 0, table[0]);
    const float y = rndQmcSobolN(uint32_t(i), 1, table[0]);
    const float z = rndQmcSobolN(uint32_t(i), 2, table[0]);
    const float w = rndQmcSobolN(uint32_t(i), 3, table[0]);

    int ix = int(x*float(AO_TILE_SIZE));
    int iy = int(y*float(AO_TILE_SIZE));
    if(ix >= AO_TILE_SIZE) ix = AO_TILE_SIZE-1;
    if(iy >= AO_TILE_SIZE) iy = AO_TILE_SIZE-1;

    uint32_t index2D = iy*AO_TILE_SIZE + ix;
    uint32_t counter = counters[index2D];
    if(counter < m_aoRaysCount)
    {
      m_aoRandomsTile[index2D*m_aoRaysCount + counter] = float2(z,w); // pack two random numbers inside target pixel storage
      counters[index2D] = counter+1;
    }
  }
}

#ifdef __ANDROID__
bool RTAO::LoadSingleMesh(const char* a_meshPath, const float* transform4x4ColMajor, AAssetManager * assetManager)
#else
bool RTAO::LoadSingleMesh(const char* a_meshPath, const float* transform4x4ColMajor)
#endif
{
  m_pAccelStruct->ClearGeom();
  
  std::cout << "[LoadScene]: mesh = " << a_meshPath << std::endl;
#if defined(__ANDROID__)
  auto currMesh = cmesh::LoadMeshFromVSGF2(assetManager, a_meshPath);
#else
  auto currMesh = cmesh::LoadMeshFromVSGF2(a_meshPath);
#endif
  if(currMesh.TrianglesNum() == 0)
  {
    std::cout << "[LoadScene]: can't load mesh '" << a_meshPath << "'" << std::endl;
    return false;
  }

  auto geomId   = m_pAccelStruct->AddGeom_Triangles3f((const float*)currMesh.vPos4f.data(), currMesh.vPos4f.size(),
                                                      currMesh.indices.data(), currMesh.indices.size(), BUILD_HIGH, sizeof(float)*4);
  float4x4 mtransform;
  if(transform4x4ColMajor != nullptr)
  {
    mtransform.set_col(0, float4(transform4x4ColMajor+0));
    mtransform.set_col(1, float4(transform4x4ColMajor+4));
    mtransform.set_col(2, float4(transform4x4ColMajor+8));
    mtransform.set_col(3, float4(transform4x4ColMajor+12));
  }

  m_pAccelStruct->ClearScene();
  m_pAccelStruct->AddInstance(geomId, mtransform);
  m_pAccelStruct->CommitScene();
  
  float3 camPos  = float3(0,0,5);
  float aspect   = float(m_width) / float(m_height);
  auto proj      = perspectiveMatrix(45.0f, aspect, 0.01f, 100.0f);
  auto worldView = lookAt(camPos, float3(0,0,0), float3(0,1,0));
  
  m_projInv      = inverse4x4(proj);
  m_worldViewInv = inverse4x4(worldView);
  m_zNearFar     = float2(0, std::numeric_limits<float>::max());
  
  m_matIdOffsets.push_back(m_matIdByPrimId.size());
  m_vertOffset.push_back(m_vNorm4f.size());
  m_matIdByPrimId.insert(m_matIdByPrimId.end(), currMesh.matIndices.begin(), currMesh.matIndices.end());
  m_triIndices.insert(m_triIndices.end(), currMesh.indices.begin(), currMesh.indices.end());
  m_vNorm4f.insert(m_vNorm4f.end(), currMesh.vNorm4f.begin(),     currMesh.vNorm4f.end());
  m_vPos4f.insert(m_vPos4f.end(), currMesh.vPos4f.begin(), currMesh.vPos4f.end());

  m_normMatrices.push_back(LiteMath::float4x4());

  //float4 boxSize = sceneBox.boxMax - sceneBox.boxMin;
  //m_aoBoxSize    = std::max(boxSize.x, std::max(boxSize.y, boxSize.z));
  m_aoBoxSize    = 100.0f;
  SetAORadius(m_presets.aoRayLength);
  m_aoRaysCount  = m_presets.aoRaysNum;

  InitQMCTable();

  return true;
}

void RTAO::SetViewport(int a_xStart, int a_yStart, int a_width, int a_height)
{
  m_width  = a_width;
  m_height = a_height;
}

void RTAO::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  auto p = timeDataByName.find(a_funcName);
  if(p == timeDataByName.end())
    return;
  a_out[0] = p->second;
}

CustomMetrics RTAO::GetMetrics() const
{
  auto traceMetrics = m_pAccelStruct->GetStats();
  CustomMetrics res = {};
  for (int i = 0; i < TREELET_ARR_SIZE; i++) {
    res.ljc_data[i] = traceMetrics.avgLJC[i];
    res.cmc_data[i] = traceMetrics.avgCMC[i];
    res.wss_data[i] = float(traceMetrics.avgWSS[i])/1000.0f;
  }
  res.common_data[0] = traceMetrics.avgNC;
  res.common_data[1] = traceMetrics.avgLC;
  res.common_data[2] = float(m_avgLCV);
  res.common_data[3] = traceMetrics.avgTC;
  res.common_data[4] = traceMetrics.avgBLB;
  res.common_data[5] = traceMetrics.avgSOC;
  res.common_data[6] = traceMetrics.avgSBL;
  res.common_data[7] = 0.0f; // traceMetrics.avgTS;
  res.size_data[0]   = traceMetrics_bvhTotalSize;
  res.size_data[1]   = traceMetrics_geomTotalSize;
  res.prims_count[0] = m_totalTrisVisiable;
  res.prims_count[1] = m_totalTris;
  return res;
}

IRenderer* MakeRTAORenderer() { return new RTAO; }
