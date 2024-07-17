#pragma once
#include "LiteMath.h"

using namespace LiteMath;

struct MaterialDataGLTF
{
  float4 baseColor;

  float metallic = 0.0f;
  float roughness = 0.0f;
  int baseColorTexId = -1;
  int metallicRoughnessTexId = -1;

  float3 emissionColor;
  int emissionTexId = -1;

  int normalTexId = -1;
  int occlusionTexId = -1;
  float alphaCutoff = 0.0;
  int alphaMode = 0;
};