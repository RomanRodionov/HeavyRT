#pragma once
#include <cstdint>
#include <memory>

#ifdef __ANDROID__
#include <android/asset_manager.h>
#endif

#include "CrossRT.h"

struct RenderPreset
{
  bool  isAORadiusInMeters;
  float aoRayLength; // in meters if isAORadiusInMeters is true else in percent of max box size
  int   aoRaysNum;
  int   numBounces;
  bool  measureOverhead;
};

#ifndef KERNEL_SLICER
struct CustomMetrics
{
  float  common_data[8];
  float  ljc_data[TREELET_ARR_SIZE];
  float  cmc_data[TREELET_ARR_SIZE];
  float  wss_data[TREELET_ARR_SIZE];   
  uint64_t size_data[2];
  uint64_t prims_count[2];
};
#endif

struct IRenderer
{
  IRenderer(){ }
  virtual ~IRenderer(){}

  virtual const char* Name() const { return ""; }

#ifdef __ANDROID__
  virtual bool LoadScene(const char* a_scenePath, AAssetManager* assetManager = nullptr) = 0;
#else
  virtual bool LoadScene(const char* a_scenePath) = 0;
#endif

#ifdef __ANDROID__
  virtual bool LoadSingleMesh(const char* a_meshPath, const float* transform4x4ColMajor, AAssetManager* assetManager = nullptr) = 0;
#else
  virtual bool LoadSingleMesh(const char* a_meshPath, const float* transform4x4ColMajor) = 0;
#endif

  virtual void Clear (uint32_t a_width, uint32_t a_height, const char* a_what) = 0;
  virtual void Render(uint32_t* imageData, uint32_t a_width, uint32_t a_height, const char* a_what, int a_passNum = 1) = 0;

  virtual void SetViewport(int a_xStart, int a_yStart, int a_width, int a_height){}
  virtual void SetPresets(const RenderPreset& a_presets){ m_presets = a_presets;}
  virtual void SetAccelStruct(std::shared_ptr<ISceneObject> a_customAccelStruct) {}
  virtual std::shared_ptr<ISceneObject> GetAccelStruct() { return nullptr; }
  
  virtual void GetExecutionTime(const char* a_funcName, float a_out[4]){}; // will be overriden in generated class

  // for future GPU impl
  //
  virtual void CommitDeviceData() {}                                     // will be overriden in generated class

  virtual void UpdateMembersPlainData() {}                               // will be overriden in generated class, optional function
  virtual void UpdateMembersVectorData() {}                              // will be overriden in generated class, optional function
  virtual void UpdateMembersTexureData() {}                              // will be overriden in generated class, optional function
  
  #ifndef KERNEL_SLICER
  virtual CustomMetrics GetMetrics() const { return CustomMetrics(); }
  #endif
  

  virtual void UpdateCamera(const LiteMath::float4x4& a_worldView, const LiteMath::float4x4& a_proj) {}

protected:

  RenderPreset m_presets;

  IRenderer(const IRenderer& rhs) {}
  IRenderer& operator=(const IRenderer& rhs) { return *this;}

  
  virtual uint32_t GetGeomNum() const  { return 0; };
  virtual uint32_t GetInstNum() const  { return 0; };
  virtual const LiteMath::float4* GetGeomBoxes() const { return nullptr; };
};

IRenderer* CreateRender(const char* a_name);
void DeleteRender(IRenderer* pImpl);

