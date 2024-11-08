#pragma once

#include <cstdint>
#include <cstddef>
#include <unordered_set>

#include "LiteMath.h"

enum BuildQuality
{
  BUILD_LOW    = 0, ///< Prefer Fast Build
  BUILD_MEDIUM = 1, ///< Standart sweep builder
  BUILD_HIGH   = 2, ///< Enable Advanced techniques like Split BVH or Early Split Clipping
  BUILD_REFIT  = 3, ///< Don't change hirarchy, recompute bouding boxes.
};

/**
\brief API to ray-scene intersection on CPU
*/
struct CRT_Hit 
{
  float    t;         ///< intersection distance from ray origin to object
  uint32_t primId; 
  uint32_t instId;
  uint32_t geomId;    ///< use 4 most significant bits for geometry type; thay are zero for triangles 
  float    coords[4]; ///< custom intersection data; for triangles coords[0] and coords[1] stores baricentric coords (u,v)
};

//#define REMAP_PRIM_ID

#define ENABLE_METRICS
constexpr int TREELET_ARR_SIZE = 5;
constexpr int treelet_sizes[TREELET_ARR_SIZE] = { 128, 512, 4096, 16384, 131072 };

#ifndef KERNEL_SLICER // CPU case only 

/**
\brief Metric statistics for traversal (average for a number of rays); 
       Note! Metrics should be evaluated only when 'ENABLE_METRICS' macro is defined!
*/
struct MetricStats
{
  float avgNC;  ///< Nodes Count
  float avgLC;  ///< Leaves Count
  float avgTC;  ///< Triangles Count
  float avgLJC[TREELET_ARR_SIZE]; ///< Long Jumps Count -- number of jumps far then treelet size during traversal;
                                  ///< please note that we also account jumps on stack.pop in this metric, howerev we can discuss this.

  float    avgCMC[TREELET_ARR_SIZE]; ///< Cache Miss Count
  uint32_t avgWSS[TREELET_ARR_SIZE]; ///< Working Set Size

  float avgBLB; ///< Bus Load in bytes
  float avgSOC; ///< Stack Operations Count (both push and pop)
  float avgSBL; ///< Stack Bytes Load (both push and pop)
  
  size_t bvhTotalSize;  ///< total byte size for all BVH data (without geometry)
  size_t geomTotalSize; ///< total byte size for all geometry data (without BVH)
};
#endif

/**
\brief API to ray-scene intersection on CPU
*/
struct ISceneObject
{
  ISceneObject(){}
  virtual ~ISceneObject(){} 
 
  /**
  \brief get implementation name  
  */
  virtual const char* Name() const = 0;

  /**
  \brief get the format name the tree build from
  */
  virtual const char* BuildName() const { return NULL; };

  /**
  \brief clear everything 
  */
  virtual void ClearGeom() = 0; 

  /**
  \brief Add geometry of type 'Triangles' to 'internal geometry library' of scene object and return geometry id
  \param a_vpos3f       - input vertex data;
  \param a_vertNumber   - vertices number. The total size of 'a_vpos4f' array is assumed to be qual to 4*a_vertNumber
  \param a_triIndices   - triangle indices (standart index buffer)
  \param a_indNumber    - number of indices, shiuld be equal to 3*triaglesNum in your mesh
  \param a_qualityLevel - bvh build quality level (low -- fast build, high -- fast traversal) 
  \param vByteStride    - byte offset from each vertex to the next one; if 0 or sizeof(float)*3 then data is tiny packed
  \return id of added geometry
  */
  virtual uint32_t AddGeom_Triangles3f(const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel = BUILD_HIGH, size_t vByteStride = sizeof(float)*3) = 0;
  
  /**
  \brief Update geometry for triangle mesh to 'internal geometry library' of scene object and return geometry id
  \param a_geomId - geometry id that should be updated. Please refer to 'AddGeom_Triangles4f' for other parameters
  
  Updates geometry. Please note that you can't: 
   * change geometry type with this fuction (from 'Triangles' to 'Spheres' for examples). 
   * increase geometry size (no 'a_vertNumber', neither 'a_indNumber') with this fuction (but it is allowed to make it smaller than original geometry size which was set by 'AddGeom_Triangles3f')
     So if you added 'Triangles' and got geom_id == 3, than you will have triangle mesh on geom_id == 3 forever and with the size you have set by 'AddGeom_Triangles3f'.
  */
  virtual void UpdateGeom_Triangles3f(uint32_t a_geomId, const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel = BUILD_HIGH, size_t vByteStride = sizeof(float)*3) = 0;
  
  /**
  \brief Clear all instances, but don't touch geometry
  */
  virtual void ClearScene() = 0; ///< 

  /**
  \brief Finish instancing and build top level acceleration structure
  */
  virtual void CommitScene(BuildQuality a_qualityLevel = BUILD_MEDIUM) = 0; ///< 
  
  /**
  \brief Add instance to scene
  \param a_geomId     - input if of geometry that is supposed to be instanced
  \param a_matrixData - float4x4 matrix, default layout is column-major

  */
  virtual uint32_t AddInstance(uint32_t a_geomId, const LiteMath::float4x4& a_matrix) = 0;
  
  /**
  \brief Add instance to scene
  \param a_instanceId
  \param a_matrixData - float4x4 matrix, the layout is column-major
  */
  virtual void     UpdateInstance(uint32_t a_instanceId, const LiteMath::float4x4& a_matrix) = 0; 
 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /**
  \brief Find nearest intersection of ray segment (Near,Far) and scene geometry
  \param posAndNear   - ray origin (x,y,z) and t_near (w)
  \param dirAndFar    - ray direction (x,y,z) and t_far (w)
  \return             - closest hit surface info
  */
  virtual CRT_Hit RayQuery_NearestHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) = 0;

  /**
  \brief Find any hit for ray segment (Near,Far). If none is found return false, else return true;
  \param posAndNear   - ray origin (x,y,z) and t_near (w)
  \param dirAndFar    - ray direction (x,y,z) and t_far (w)
  \return             - true if a hit is found, false otherwaise
  */
  virtual bool    RayQuery_AnyHit(LiteMath::float4 posAndNear, LiteMath::float4 dirAndFar) = 0;
  
  #ifndef KERNEL_SLICER
  virtual MetricStats GetStats();
  virtual void ResetStats();
  #endif

  #ifdef ENABLE_METRICS
  void     ResetVarLC() { m_stats.LC2 = 0; } 
  uint64_t GetVarLC  () const { return m_stats.LC2; }
  #endif

  virtual uint32_t GetGeomNum() const  { return 0; };
  virtual uint32_t GetInstNum() const  { return 0; };
  virtual const LiteMath::float4* GetGeomBoxes() const { return nullptr; };

protected:
  
  #ifndef KERNEL_SLICER
  struct MetricsI
  {
    uint64_t NC;  ///< Nodes Count
    uint64_t LC;  ///< Leaves Count
    uint64_t LC2; ///< Leaves Count which we can reset and retrieve per ray 
    uint64_t TC;  ///< Triangles Count
    uint64_t LJC[TREELET_ARR_SIZE]; ///< Long Jumps Count -- number of jumps far then treelet size during traversal;
                                    ///< please note that we also account jumps on stack.pop in this metric, howerev we can discuss this.

    uint32_t                     CMC[TREELET_ARR_SIZE]; // Cache Miss Count
    std::unordered_set<uint64_t> WSS[TREELET_ARR_SIZE]; // Working Set Size

    uint64_t BLB; ///< Bus Load in bytes
    uint64_t SOC; ///< Stack Operations Count (both push and pop)
    uint64_t SBL; ///< Stack Bytes Load (both push and pop)
    uint32_t raysNumber;  ///< already traversed rays number 
    
    size_t bvhTotalSize;  ///< total byte size for all BVH data (without geometry)
    size_t geomTotalSize; ///< total byte size for all geometry data (without BVH)

    void clear() { MetricsI zerows = {}; *this = zerows; }
  } m_stats;
  #endif
};

ISceneObject* CreateSceneRT(const char* a_implName, const char* a_buildName, const char* a_layoutName);
void DeleteSceneRT(ISceneObject* a_impl);
