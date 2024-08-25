#include "rtao.h"
#include "render_common.h"

#include <cfloat>
#include <chrono>

using namespace LiteMath;

static inline float attenuationFunc(float r, float R) {
//  return std::clamp(r/R, 0.0f, 1.0f); // HBAO attenuation function
  return 0.0f;
}

static inline uint32_t float4ToRGBA(LiteMath::float4 v)
{
  LiteMath::int4 intV = LiteMath::int4(clamp(v, 0, 1)*255);
  uint32_t rgba = intV.x + (intV.y<<8) + (intV.z<<16) + (intV.w<<24);
  return rgba;
}

static inline uint32_t floatToGrayRGBA(float v)
{
  return float4ToRGBA(LiteMath::float4(v,v,v, 1.0f));
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

void RTAO::CalcAOBlock(uint32_t* a_outColor, uint32_t a_width, uint32_t a_height, uint32_t a_passNumber)
{
  auto start = std::chrono::high_resolution_clock::now();
  
  #ifndef _DEBUG
  #ifndef ENABLE_METRICS
  #pragma omp parallel for collapse (2)
  #endif
  #endif
  for (int j = 0; j < int(a_height); ++j)
    for (int i = 0; i < int(a_width); ++i)
      for(uint32_t k=0;k<a_passNumber;k++)
        CalcAO(a_outColor, i, j);

  timeDataByName["CalcAOBlock"] = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.f;
}

void RTAO::CalcAO(uint32_t* a_outColor, uint32_t tidX, uint32_t tidY)
{
  float4 rayPosAndNear, rayDirAndFar, hitPosNorm;
  float visibility;

  kernel_InitEyeRay  (tidX, tidY, &rayPosAndNear, &rayDirAndFar, &visibility); // ==> (rayPosAndNear, rayDirAndFar, visibility)
  kernel_TraceEyeRay2(tidX, tidY, &rayPosAndNear, &rayDirAndFar, &hitPosNorm); // ==> (hitPos,hitNorm)
  
  for(uint32_t tidZ = 0; tidZ < m_aoRaysCount; tidZ++) {
    kernel_InitAORay (tidX, tidY, tidZ, &rayPosAndNear, &rayDirAndFar, &hitPosNorm); // ==> (rayPosAndNear, rayDirAndFar, cosAlpha)
    kernel_TraceAORay(tidX, tidY,       &rayPosAndNear, &rayDirAndFar, &visibility);        // ==> visibility
  }

  kernel_AO2Color(tidX, tidY, &hitPosNorm, &visibility, a_outColor); // ==> a_outColor
}

void RTAO::kernel_AO2Color(uint32_t tidX, uint32_t tidY, const float4* positions, const float* visibility, uint32_t* out_color)
{
  if(positions->y >= AO_HIT_BACK)
    out_color[tidY*m_width + tidX] = 0;
  else
  {
    const float sampleCount = float(AO_PASS_COUNT * m_aoRaysCount);
    const float normCoeff   = 1.0f / sampleCount;
    const float resAO       = std::pow((*visibility) * normCoeff, m_power);
    out_color[tidY*m_width + tidX] = floatToGrayRGBA(resAO);
  }
}

static float3 calcTriangleNormal(const float3 *A_pos, const float3 *B_pos, const float3 *C_pos) {
  //CounterClockwise normal
  return cross((*B_pos) - (*A_pos), (*C_pos) - (*A_pos));
}

void RTAO::kernel_TraceEyeRay2(uint32_t tidX, uint32_t tidY, const float4* rayPosAndNear,
                               const float4* rayDirAndFar, float4* positions)
{
  const float4 rayPos = *rayPosAndNear;
  const float4 rayDir = *rayDirAndFar ;

  CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(rayPos, rayDir);

  if(hit.geomId != uint32_t(-1))
  {
    const float2 uv       = float2(hit.coords[0], hit.coords[1]);
    const float3 hitPos   = to_float3(rayPos) + 0.9999995f*hit.t*to_float3(rayDir);

    const uint triOffset  = m_matIdOffsets[hit.geomId];
    const uint vertOffset = m_vertOffset  [hit.geomId];

    const uint A = m_triIndices[(triOffset + hit.primId)*3 + 0];
    const uint B = m_triIndices[(triOffset + hit.primId)*3 + 1];
    const uint C = m_triIndices[(triOffset + hit.primId)*3 + 2];

    const float3 A_pos = to_float3(m_vPos4f[A + vertOffset]);
    const float3 B_pos = to_float3(m_vPos4f[B + vertOffset]);
    const float3 C_pos = to_float3(m_vPos4f[C + vertOffset]);

    float3 triangleNormal = calcTriangleNormal(&A_pos, &B_pos, &C_pos);

    const float3 A_norm = to_float3(m_vNorm4f[A + vertOffset]);
    const float3 B_norm = to_float3(m_vNorm4f[B + vertOffset]);
    const float3 C_norm = to_float3(m_vNorm4f[C + vertOffset]);

    float3 hitNorm      = (1.0f - uv.x - uv.y)*A_norm + uv.y*B_norm + uv.x*C_norm;
   
    // transform surface point with matrix and flip normal if needed
    //
    hitNorm              = normalize(mymul3x3(m_normMatrices[hit.instId], hitNorm));
    triangleNormal       = normalize(mymul3x3(m_normMatrices[hit.instId], triangleNormal));

    if(dot(triangleNormal,hitNorm) < 0.0f)
      hitNorm *= -1.0f;

    const float dotp     = dot(to_float3(rayDir), triangleNormal);
    const float flipNorm = (dotp > 1e-5f) ? -1.0f : 1.0f; // beware of transparent materials which use normal sign to identity "inside/outside" glass for example
    const uint normalCompressed = encodeNormal(hitNorm*flipNorm);

    *positions = to_float4(hitPos, as_float(normalCompressed));
  }
  else
  {
    const uint normalCompressed = encodeNormal(float3(0,1,0));
    *positions = float4(0.0f, AO_HIT_BACK, 0.0f, as_float(normalCompressed));
  }
}

void RTAO::kernel_InitEyeRay(uint32_t tidX, uint32_t tidY, float4* rayPosAndNear, float4* rayDirAndFar, float* visibility)
{
  float3 rayDir = EyeRayDirNormalized((float(tidX)+0.5f)/float(m_width), 
                                      (float(tidY)+0.5f)/float(m_height), m_projInv);
  float3 rayPos = float3(0,0,0);

  transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
  
  *rayPosAndNear = to_float4(rayPos, m_zNearFar.x); // 0.0f
  *rayDirAndFar  = to_float4(rayDir, m_zNearFar.y); // FLT_MAX
  *visibility    = 0.0f;
}


void RTAO::kernel_InitAORay(uint32_t tidX, uint32_t tidY, uint32_t tidZ,
                            float4* rayPosAndNear, float4* rayDirAndFar, const float4* positions)
{
  const float4 rayPos = *positions;
  if (rayPos.y>=AO_HIT_BACK) // no hit point of screen
    return;

  const float3 normal = decodeNormal(as_uint(rayPos.w)); // to_float3(*normals);

  const uint32_t xTiled = tidX % AO_TILE_SIZE;
  const uint32_t yTiled = tidY % AO_TILE_SIZE;

  float2 uv      = m_aoRandomsTile[(yTiled * AO_TILE_SIZE + xTiled)*m_aoRaysCount + tidZ];
  float3 rayDir  = MapSampleToCosineDistribution(uv.x, uv.y, normal, normal, 1.0f);
  float3 rayPos2 = OffsRayPos(to_float3(rayPos), normal, rayDir);

  *rayPosAndNear = to_float4(rayPos2, 0.0f); // rayPos.w
  *rayDirAndFar  = to_float4(rayDir, m_aoMaxRadius);
}

void RTAO::kernel_TraceAORay(uint32_t tidX, uint32_t tidY,
                             const float4* rayPosAndNear, const float4* rayDirAndFar, float* out_visibility)
{
  const float4 rayPos = *rayPosAndNear;
  const float4 rayDir = *rayDirAndFar ;
  if (rayPos.y>=AO_HIT_BACK)
    return;
  
  #ifdef RTAO_USE_CLOSEST_HIT
  CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(rayPos, rayDir);
  float visibility = 1.0f;
  if (hit.geomId != uint32_t(-1)) {
    visibility = attenuationFunc(hit.t, m_aoMaxRadius);
  }
  *out_visibility = *out_visibility + visibility;
  #else
  bool hit = m_pAccelStruct->RayQuery_AnyHit(rayPos, rayDir);
  if(hit)
    *out_visibility = *out_visibility + 1.0f;
  #endif
}
