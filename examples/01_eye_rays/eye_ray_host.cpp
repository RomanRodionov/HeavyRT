#include <cfloat>
#include <cstring>
#include <sstream>

#include "eye_ray.h"
#include "loader/hydraxml.h"
#include "loader/cmesh.h"
#include "loader/gltf_loader.h"
#include "gltf_data.h"
#include "Timer.h"

using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;

using LiteMath::perspectiveMatrix;
using LiteMath::lookAt;
using LiteMath::inverse4x4;

EyeRayCaster::EyeRayCaster() 
{ 
  m_pAccelStruct = nullptr;
}

void EyeRayCaster::SetViewport(int a_xStart, int a_yStart, int a_width, int a_height)
{
  m_width  = a_width;
  m_height = a_height;
  m_packedXY.resize(m_width*m_height);
}

#if defined(__ANDROID__)
bool EyeRayCaster::LoadScene(const char* a_scenePath, AAssetManager* assetManager)
#else
bool EyeRayCaster::LoadScene(const char* a_scenePath)
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
bool EyeRayCaster::LoadSceneHydra(const std::string& a_path, AAssetManager* assetManager)
#else
bool EyeRayCaster::LoadSceneHydra(const std::string& a_path)
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

    m_projInv      = inverse4x4(proj);
    m_worldViewInv = inverse4x4(worldView);

    break; // take first cam
  }

  std::vector<uint64_t> trisPerObject;
  trisPerObject.reserve(1000);
  m_totalTris = 0;
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
  }
  
  m_totalTrisVisiable = 0;
  m_pAccelStruct->ClearScene();
  for(auto inst : scene.InstancesGeom())
  {
    m_pAccelStruct->AddInstance(inst.geomId, inst.matrix);
    m_totalTrisVisiable += trisPerObject[inst.geomId];
  }
  m_pAccelStruct->CommitScene();

  return true;
}


struct DataRefs
{
  DataRefs(std::unordered_map<int, std::pair<uint32_t, BBox3f>> &a_loadedMeshesToMeshId, std::vector<uint64_t>& a_trisPerObject, Box4f& a_sceneBBox,
           int& a_gltfCam, float4x4& a_worldViewInv, std::shared_ptr<ISceneObject> a_pAccelStruct, uint64_t& a_totalTris) : 
           m_loadedMeshesToMeshId(a_loadedMeshesToMeshId), m_trisPerObject(a_trisPerObject), m_sceneBBox(a_sceneBBox),
           m_gltfCamId(a_gltfCam), m_worldViewInv(a_worldViewInv), m_pAccelStruct(a_pAccelStruct), m_totalTris(a_totalTris) {}

  
  std::unordered_map<int, std::pair<uint32_t, BBox3f>>& m_loadedMeshesToMeshId;
  std::vector<uint64_t>&                                m_trisPerObject;
  Box4f&                                                m_sceneBBox;

  int&                          m_gltfCamId;
  float4x4&                     m_worldViewInv;
  std::shared_ptr<ISceneObject> m_pAccelStruct;
  uint64_t&                     m_totalTris;
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
      }
    }
    auto tmp_box = a_refs.m_loadedMeshesToMeshId[a_node.mesh].second;
    Box4f mesh_box {};
    mesh_box.boxMin = to_float4(tmp_box.boxMin, 1.0);
    mesh_box.boxMax = to_float4(tmp_box.boxMax, 1.0);
    Box4f inst_box = {};
    inst_box.boxMin = nodeMatrix * mesh_box.boxMin;
    inst_box.boxMax = nodeMatrix * mesh_box.boxMax;

    a_refs.m_sceneBBox.include(inst_box);
    a_refs.m_pAccelStruct->AddInstance(a_refs.m_loadedMeshesToMeshId[a_node.mesh].first, nodeMatrix);
  }
}


#ifdef __ANDROID__
bool EyeRayCaster::LoadSceneGLTF(const std::string& a_path, AAssetManager* assetManager)
#else
bool EyeRayCaster::LoadSceneGLTF(const std::string& a_path)
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
  for(size_t i = 0; i < gltfModel.cameras.size(); ++i)
  {
    auto gltfCam = gltfModel.cameras[i];
    if(gltfCam.type == "perspective")
    {
      auto proj = perspectiveMatrix(gltfCam.perspective.yfov / DEG_TO_RAD, aspect,
                                    gltfCam.perspective.znear, gltfCam.perspective.zfar);
      m_projInv = inverse4x4(proj);

      m_gltfCamId = i;
      break; // take first perspective cam
    }
    else if (gltfCam.type == "orthographic")
    {
      std::cerr << "Orthographic camera not supported!" << std::endl;
    }
  }


  // load and instance geometry
  std::vector<uint64_t> trisPerObject;
  trisPerObject.reserve(1000);

  m_totalTris = 0;
  m_pAccelStruct->ClearGeom();
  m_pAccelStruct->ClearScene();

  Box4f sceneBBox;
  std::unordered_map<int, std::pair<uint32_t, BBox3f>> loaded_meshes_to_meshId;
  for(size_t i = 0; i < scene.nodes.size(); ++i)
  {
    const tinygltf::Node node = gltfModel.nodes[scene.nodes[i]];
    auto identity = LiteMath::float4x4();
    LoadGLTFNodesRecursive(gltfModel, node, identity, 
                           DataRefs(loaded_meshes_to_meshId, trisPerObject, sceneBBox, m_gltfCamId, m_worldViewInv, m_pAccelStruct, m_totalTris));
  }

  // glTF scene can have no cameras specified
  if(m_gltfCamId == -1)
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

  return true;
}


#if defined(__ANDROID__)
bool EyeRayCaster::LoadSingleMesh(const char* a_meshPath, const float* transform4x4ColMajor, AAssetManager * assetManager)
#else
bool EyeRayCaster::LoadSingleMesh(const char* a_meshPath, const float* transform4x4ColMajor)
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

  return true;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

void EyeRayCaster::Render(uint32_t* a_outColor, uint32_t a_width, uint32_t a_height, const char* a_what, int a_passNum)
{
  CastRaySingleBlock(a_width*a_height, a_outColor, a_passNum);
}

void EyeRayCaster::CastRaySingleBlock(uint32_t tidX, uint32_t * out_color, uint32_t a_numPasses)
{
  profiling::Timer timer;
  
  #ifndef _DEBUG
  #ifndef ENABLE_METRICS
  #pragma omp parallel for default(shared)
  #endif
  #endif
  for(int i=0;i<tidX;i++)
    CastRaySingle(i, out_color);

  //printf("CastSingleRayBlock...\n");
  //
  //auto a_width  = tidX;
  //auto a_height = tidY;
  //m_avgLCV = 0.0;
  //
  //#ifndef _DEBUG
  //#ifndef ENABLE_METRICS
  //#pragma omp parallel for collapse (2)
  //#endif
  //#endif
  //for (int j = 0; j < int(a_height); j+=BSIZE)
  //  for (int i = 0; i < int(a_width); i+=BSIZE)
  //    CastRayPacket(i, j, out_color);
  //
  //m_avgLCV /= double( (a_height*a_width)/(BSIZE*BSIZE) );
  //
  timeDataByName["CastRaySingleBlock"] = timer.getElapsedTime().asMilliseconds();
  //printf("CastRaySingleBlock %8.3f sec\n", timer.getElapsedTime().asSeconds());
}

const char* EyeRayCaster::Name() const
{
  std::stringstream strout;
  strout << "EyeRayCaster(" << m_pAccelStruct->Name() << ")";
  m_tempName = strout.str();
  return m_tempName.c_str();
}

CustomMetrics EyeRayCaster::GetMetrics() const 
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
  res.size_data[0] = traceMetrics.bvhTotalSize;
  res.size_data[1] = traceMetrics.geomTotalSize;
  res.prims_count[0] = m_totalTrisVisiable;
  res.prims_count[1] = m_totalTris;
  return res;
}


void EyeRayCaster::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  auto p = timeDataByName.find(a_funcName);
  if(p == timeDataByName.end())
    return;
  a_out[0] = p->second;
}


IRenderer* MakeEyeRayShooterRenderer(const char* a_name) { return new EyeRayCaster; }

