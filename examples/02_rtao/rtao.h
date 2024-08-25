#pragma once

#include <cstdint>
#include <memory>
#include <array>
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>

#include "LiteMath.h"
#include "CrossRT.h"
#include "IRenderer.h"

#define RTAO_USE_CLOSEST_HIT

class RTAO : public IRenderer
{
public:
  RTAO();

  virtual void SceneRestrictions(uint32_t a_restrictions[4]) const
  {
    uint32_t maxMeshes            = 1024;
    uint32_t maxTotalVertices     = 8'000'000;
    uint32_t maxTotalPrimitives   = 8'000'000;
    uint32_t maxPrimitivesPerMesh = 4'000'000;

    a_restrictions[0] = maxMeshes;
    a_restrictions[1] = maxTotalVertices;
    a_restrictions[2] = maxTotalPrimitives;
    a_restrictions[3] = maxPrimitivesPerMesh;
  }

  virtual const char* Name() const override { return "RTAO"; }

#ifdef __ANDROID__
  virtual bool LoadScene(const char* a_scenePath, AAssetManager* assetManager = nullptr) override;
#else
  virtual bool LoadScene(const char* a_scenePath) override;
#endif

#ifdef __ANDROID__
  bool LoadSingleMesh(const char* a_meshPath, const float* transform4x4ColMajor, AAssetManager* assetManager = nullptr) override;
#else
  bool LoadSingleMesh(const char* a_meshPath, const float* transform4x4ColMajor) override;
#endif

  void Clear (uint32_t a_width, uint32_t a_height, const char* a_what) override { }
  void Render(uint32_t* imageData, uint32_t a_width, uint32_t a_height, const char* a_what, int a_passNum) override;
  void SetViewport(int a_xStart, int a_yStart, int a_width, int a_height) override; //ignores start arguments
  void SetAccelStruct(std::shared_ptr<ISceneObject> a_customAccelStruct) override { m_pAccelStruct = a_customAccelStruct;}
  std::shared_ptr<ISceneObject> GetAccelStruct() override { return m_pAccelStruct; }
  
  void GetExecutionTime(const char* a_funcName, float a_out[4]) override;
  void CommitDeviceData() override {}

  #ifndef KERNEL_SLICER
  CustomMetrics GetMetrics() const override;
  #endif

  void UpdateCamera(const LiteMath::float4x4& a_worldView, const LiteMath::float4x4& a_proj) override 
  {
    m_projInv      = inverse4x4(a_proj);
    m_worldViewInv = inverse4x4(a_worldView);
  }

  virtual uint32_t GetGeomNum() const  { return m_pAccelStruct->GetGeomNum(); };
  virtual uint32_t GetInstNum() const  { return m_pAccelStruct->GetInstNum(); };
  virtual const LiteMath::float4* GetGeomBoxes() const { return m_pAccelStruct->GetGeomBoxes(); };

protected:

  void SetAORadius(float radius);

  void kernel_InitEyeRay(uint32_t tidX, uint32_t tidY, LiteMath::float4* rayPosAndNear, LiteMath::float4* rayDirAndFar, float* visibility);
  void kernel_AO2Color(uint32_t tidX, uint32_t tidY, const LiteMath::float4* positions, const float* visibility, uint32_t* out_color);

  void kernel_TraceEyeRay2(uint32_t tidX, uint32_t tidY, const LiteMath::float4* rayPosAndNear,
                           const LiteMath::float4* rayDirAndFar, LiteMath::float4* positions);

  #ifdef KERNEL_SLICER
  virtual void CalcAO(uint32_t* a_outColor __attribute__((size("tidX*tidY"))), uint32_t tidX, uint32_t tidY);
  #else
  virtual void CalcAO(uint32_t* a_outColor, uint32_t tidX, uint32_t tidY);
  #endif

  virtual void CalcAOBlock(uint32_t* a_outColor, uint32_t a_width, uint32_t a_height, uint32_t a_passNumber);

  void CastSingleAORay(uint32_t tidX, uint32_t tidY, uint32_t tidZ, float* ao_tex, int a_passNumber);
  void kernel_InitAORay(uint32_t tidX, uint32_t tidY, uint32_t tidZ,
                        LiteMath::float4* rayPosAndNear, LiteMath::float4* rayDirAndFar,
                        const LiteMath::float4* positions);
  void kernel_TraceAORay(uint32_t tidX, uint32_t tidY,
                         const LiteMath::float4* rayPosAndNear, const LiteMath::float4* rayDirAndFar, float* out_visibility);

  static constexpr uint32_t AO_PASS_COUNT = 1;
  static constexpr uint32_t AO_TILE_SIZE  = 128;
  static constexpr float    AO_HIT_BACK   = 10000000.0f;

  void InitQMCTable();

#ifdef __ANDROID__
  bool LoadSceneHydra(const std::string& a_path, AAssetManager* assetManager = nullptr);
  bool LoadSceneGLTF(const std::string& a_path, AAssetManager* assetManager = nullptr);
#else
  bool LoadSceneHydra(const std::string& a_path);
  bool LoadSceneGLTF(const std::string& a_path);
#endif

  uint32_t m_width;
  uint32_t m_height;
  float    m_aoMaxRadius = 2.0f;
  float    m_power       = 1.0f;
  float    m_aoBoxSize   = 100.0f; 
  uint32_t m_aoRaysCount = 4;

  uint64_t m_totalTris         = 0;
  uint64_t m_totalTrisVisiable = 0;
  double   m_avgLCV = 0.0;
  uint64_t traceMetrics_bvhTotalSize  = 0;
  uint64_t traceMetrics_geomTotalSize = 0;

  LiteMath::float2   m_zNearFar;
  LiteMath::float4x4 m_projInv;
  LiteMath::float4x4 m_worldViewInv;
  std::vector<LiteMath::float4x4> m_normMatrices; ///< per instance normal matrix, local to world

  std::shared_ptr<ISceneObject> m_pAccelStruct;

  std::vector<uint32_t>        m_matIdOffsets;  ///< offset = m_matIdOffsets[geomId]
  std::vector<uint32_t>        m_matIdByPrimId; ///< matId  = m_matIdByPrimId[offset + primId]
  std::vector<uint32_t>        m_triIndices;    ///< (A,B,C) = m_triIndices[(offset + primId)*3 + 0/1/2]

  std::vector<uint32_t>          m_vertOffset;  ///< vertOffs = m_vertOffset[geomId]
  std::vector<LiteMath::float4>  m_vNorm4f;     ///< vertNorm = m_vNorm4f[vertOffs + vertId]
  std::vector<LiteMath::float4>  m_vPos4f;
  
  std::vector<LiteMath::float2>  m_aoRandomsTile;

  // color palette to select color for objects based on mesh/instance id
  static constexpr uint32_t palette_size = 20;
  static constexpr uint32_t m_palette[palette_size] = {
    0xffe6194b, 0xff3cb44b, 0xffffe119, 0xff0082c8,
    0xfff58231, 0xff911eb4, 0xff46f0f0, 0xfff032e6,
    0xffd2f53c, 0xfffabebe, 0xff008080, 0xffe6beff,
    0xffaa6e28, 0xfffffac8, 0xff800000, 0xffaaffc3,
    0xff808000, 0xffffd8b1, 0xff000080, 0xff808080
  };

  std::unordered_map<std::string, float> timeDataByName;
};


