#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <chrono>  // C++11
//#if 1 // Ose
//#include <omp.h>
//#endif
#include "CrossRT.h"
#include "raytrace_common.h"
#include "nanort/nanort.h"
#include "NanoRText.h"

using LiteMath::dot;
using LiteMath::sign;
using LiteMath::cross;
using LiteMath::float4x4;
using LiteMath::uint2;
using LiteMath::int2;
using LiteMath::float2;
using LiteMath::float3;
using LiteMath::float4;
using LiteMath::normalize;
using LiteMath::inverse4x4;
using LiteMath::to_float3;
using LiteMath::Box4f;


struct NanoRT : public ISceneObject
{
  NanoRT(){}
  ~NanoRT() override {}

  const char* Name() const override { return "NanoRT"; }

  void ClearGeom() override;

  uint32_t AddGeom_Triangles3f(const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride) override;
  void     UpdateGeom_Triangles3f(uint32_t a_geomId, const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride) override;

  void ClearScene() override;
  void CommitScene(BuildQuality a_qualityLevel) override;

  uint32_t AddInstance(uint32_t a_geomId, const float4x4& a_matrix) override;
  void     UpdateInstance(uint32_t a_instanceId, const float4x4& a_matrix) override;

  CRT_Hit  RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar) override;
  bool     RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar) override;

  uint32_t GetGeomNum() const override { return uint32_t(m_geomBoxes.size()); }
  uint32_t GetInstNum() const override { return uint32_t(m_instBoxes.size()); }
  const LiteMath::float4* GetGeomBoxes() const override { return (const LiteMath::float4*)m_geomBoxes.data(); }

protected:

  std::vector< nanort::BVHAccel<float> > m_accels;

  std::vector<Box4f> m_geomBoxes;
  std::vector<Box4f> m_instBoxes;

  std::vector<float4x4> m_instMatricesInv; ///< inverse instance matrices
  std::vector<float4x4> m_instMatricesFwd; ///< instance matrices

  std::vector<float>    m_verts;
  std::vector<float4>   m_vertPos;
  std::vector<uint32_t> m_indices;

  std::vector<uint2>    m_geomOffsets;
  std::vector<uint32_t> m_geomIdByInstId;

  std::vector< std::vector<BoxHit> > m_wksp0_thr;
};


struct NanoRTExt : public NanoRT
  {
  ~NanoRTExt() override {}

  const char* Name() const override { return "NanoRTExt"; }

  void CommitScene(BuildQuality a_qualityLevel) override;

  CRT_Hit  RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar) override;

#if 1 // ose
  /// Provides tree optimization (re-arrangement in memory) or not. No by default
  virtual bool CanDoTreeOptimization() const { return true; };
  /// Prints data for optimization of BVH tree, if available
  virtual int PrintTreeOptimizData(const char *fname) const;
  /// Optimize tree, if has data for it
  virtual int OptimizeBVH();

  /// Get statistics
  virtual MetricStats GetStats();


  protected:
    /// Determines which nodes (of a sub-tree) constitute its root cluster
    /// and which (those "below" the cluster in the tree hierarchy) will be
    /// the roots of child clusters.
    int RootCluster4Subtree(int root_node_idx,
      const std::vector<double>& w0,
      double threshold,
      int max_num_nodes,
      std::vector<int>& rcluster_idx,
      std::vector<int>& child_cluster_root_idx) const;

    /// Suggests ordering of the nodes of a subtree "by clusters". First comes the
    /// "root cluster" then "child clusters". Each of the latter, in turn, consists
    /// of its root cluster, past which come child clusters, of which... etc.
    int ClusterizeSubtree(int root_node_idx,
      const std::vector<double>& w0,
      double threshold,
      std::vector<int>& indices) const;

  public:
    /// Changes the order of the tree nodes in memory WHLE RETAINING THE TREE TOPOLOGY
    /// and updates the links to the children kept in the nodes accordingly.
    int ReArrangeByClusters(const std::vector<double>& w0, double threshold);

    /// Traverse subtree and for small-size nodes calculate frequency of hits by scaling that of parent
    int ImproveSmallNodesProb(int root_node_idx, std::vector< real3<double> >& num_hits) const;
#endif

  protected:
    /// Number of hits by node walls. Outer index is that of scene geometry instance (=of BVHAccelExt object), then index of node, then index of wall
    std::vector < std::vector< real3<double> > > m_num_hits;

    /// Weight of node. Outer index is that of scene geometry instance (=of BVHAccelExt object)
    std::vector < std::vector<double> > m_w0;
    /// Workspace e.g. for calculation of "weights" of nodes in ray tracing. The outer index is thread ID
    std::vector < std::vector<int> > m_wksp1_thr;
    /// Workspace e.g. for calculation of "weights" of nodes in ray tracing. The outer index is thread ID
    std::vector < std::vector<int> > m_wksp2_thr;
  };


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr size_t reserveSize = 1000;

void NanoRT::ClearGeom()
{
  m_accels.reserve(std::max(reserveSize, m_accels.capacity()));
  m_accels.resize(0);

  m_verts.reserve(std::max<size_t>(100000 * 3, m_verts.capacity()));
  m_vertPos.reserve(std::max<size_t>(100000, m_vertPos.capacity()));
  m_indices.reserve(std::max<size_t>(100000 * 3, m_indices.capacity()));

  m_verts.resize(0);
  m_vertPos.resize(0);
  m_indices.resize(0);

  m_geomOffsets.reserve(std::max(reserveSize, m_geomOffsets.capacity()));
  m_geomOffsets.resize(0);

  m_geomBoxes.reserve(std::max<size_t>(reserveSize, m_geomBoxes.capacity()));
  m_geomBoxes.resize(0);

  ClearScene();
}

uint32_t NanoRT::AddGeom_Triangles3f(const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride)
{
  const size_t vStride = vByteStride / 4;
  assert(vByteStride % 4 == 0);

  const uint32_t currGeomId = uint32_t(m_geomOffsets.size());
  const size_t oldSizeVert  = m_vertPos.size();
  const size_t oldSizeInd   = m_indices.size();

  m_geomOffsets.push_back(uint2(oldSizeInd, oldSizeVert));

  m_vertPos.resize(oldSizeVert + a_vertNumber);
  m_indices.resize(oldSizeInd + a_indNumber);

  Box4f bbox;
  for (size_t i = 0; i < a_vertNumber; i++)
  {
    const float4 v = float4(a_vpos3f[i * vStride + 0], a_vpos3f[i * vStride + 1], a_vpos3f[i * vStride + 2], 1.0f);
    m_vertPos[oldSizeVert + i] = v;
    bbox.include(v);
  }

  m_geomBoxes.push_back(bbox);

  for(size_t i = 0 ; i < a_indNumber; i++)
    m_indices[oldSizeInd + i] = a_triIndices[i];

  const size_t old_verts_size  = m_verts.size();
  m_verts.resize(old_verts_size + a_vertNumber * 3);
  for(size_t y = 0; y < a_vertNumber; y++) 
  {
    m_verts[old_verts_size + (y * 3) + 0] = a_vpos3f[y * 4 + 0];
    m_verts[old_verts_size + (y * 3) + 1] = a_vpos3f[y * 4 + 1];
    m_verts[old_verts_size + (y * 3) + 2] = a_vpos3f[y * 4 + 2];
  }

  std::cout << "[NanoRT: Build BVH] " << std::endl;

  nanort::BVHBuildOptions<float> build_options;  // Use default option
  build_options.cache_bbox = false;

  auto t_start = std::chrono::system_clock::now();

  nanort::TriangleMesh<float> triangle_mesh(m_verts.data() + old_verts_size, m_indices.data() + oldSizeInd, sizeof(float) * 3);
  nanort::TriangleSAHPred<float> triangle_pred(m_verts.data() + old_verts_size, m_indices.data() + oldSizeInd, sizeof(float) * 3);

  const size_t n_trgs = a_indNumber / 3;
  printf("num_triangles = %lu\n", n_trgs);

  nanort::BVHAccel<float> accel;
  bool ret = accel.Build(n_trgs, triangle_mesh, triangle_pred, build_options);
  if(!ret)
  {
    std::cout << "nanort::accel.Build is failed" << std::endl;
    exit(0);
  }

  auto t_end = std::chrono::system_clock::now();

  std::chrono::duration<double, std::milli> ms = t_end - t_start;
  std::cout << "BVH build time: " << ms.count() << " [ms]\n";

  nanort::BVHBuildStatistics stats = accel.GetStatistics();

  printf("  BVH statistics:\n");
  printf("    # of leaf   nodes: %d\n", stats.num_leaf_nodes);
  printf("    # of branch nodes: %d\n", stats.num_branch_nodes);
  printf("  Max tree depth     : %d\n", stats.max_tree_depth);
  float bmin[3], bmax[3];
  accel.BoundingBox(bmin, bmax);
  printf("  Bmin               : %f, %f, %f\n", bmin[0], bmin[1], bmin[2]);
  printf("  Bmax               : %f, %f, %f\n", bmax[0], bmax[1], bmax[2]);

  m_accels.push_back(accel);
  return currGeomId;
}

void NanoRT::UpdateGeom_Triangles3f(uint32_t a_geomId, const float* a_vpos3f, size_t a_vertNumber, const uint32_t* a_triIndices, size_t a_indNumber, BuildQuality a_qualityLevel, size_t vByteStride)
{
  std::cout << "[NanoRT::UpdateGeom_Triangles3f]: " 
            << "not implemeted!" << std::endl;
}

void NanoRT::ClearScene()
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

void NanoRT::CommitScene(BuildQuality a_qualityLevel)
{
  // Initialize workspace. 
  int num_threads = 1;

  ///#ifndef _DEBUG
  //num_threads = omp_get_max_threads();
  ///#endif

  m_wksp0_thr.reserve(num_threads);
  m_wksp0_thr.resize(num_threads);

  // reset stats
  //  
  m_stats.clear();
}

uint32_t NanoRT::AddInstance(uint32_t a_geomId, const float4x4& a_matrix)
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

void NanoRT::UpdateInstance(uint32_t a_instanceId, const float4x4& a_matrix)
{
  std::cout << "[NanoRT::UpdateInstance]: " 
            << "not implemeted!" << std::endl;
}

CRT_Hit NanoRT::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  int thread_id = 0;
  ///#ifndef _DEBUG
  //thread_id = omp_get_thread_num();
  ///#endif

  #ifdef ENABLE_METRICS
  m_stats.raysNumber++;
  #endif

  float3 invRayDir = SafeInverse(to_float3(dirAndFar));
  float3 rayPos    = to_float3(posAndNear);
  
  std::vector<BoxHit>& boxMinHits = m_wksp0_thr[thread_id];
  boxMinHits.clear();
  boxMinHits.reserve(64);
  
  const float tMin = posAndNear.w;
  const float tMax = dirAndFar.w; 

  // (1) check all instance boxes and all triangles inside box. 
  for (uint32_t i = 0; i < uint32_t(m_instBoxes.size()); i++)
  {
    const float2 minMax = RayBoxIntersection2(rayPos, invRayDir, to_float3(m_instBoxes[i].boxMin), to_float3(m_instBoxes[i].boxMax));
    if( (minMax.x <= minMax.y) && (minMax.y >= tMin) && (minMax.x <= tMax) )
      boxMinHits.push_back(make_BoxHit(i, minMax.x));
  }
  
  // (2) sort all hits by hit distance to process nearest first
  //
  std::sort(boxMinHits.begin(), boxMinHits.end(), [](const BoxHit a, const BoxHit b){ return a.tHit < b.tHit; });

  // (3) process all potential hits
  //
  CRT_Hit hit;
  hit.t      = tMax;
  hit.primId = uint32_t(-1);
  hit.instId = uint32_t(-1);
  hit.geomId = uint32_t(-1);

  for (uint32_t boxId = 0; boxId < uint32_t(boxMinHits.size()); boxId++)
  {
    if (boxMinHits[boxId].tHit > hit.t) // already found hit that is closer than bounding box hit
      break;

    const uint32_t instId = boxMinHits[boxId].id;
    const uint32_t geomId = m_geomIdByInstId[instId];
    const uint2 startSize = m_geomOffsets[geomId];

    // transform ray with matrix to local space
    //
    float3 ray_pos = mul4x3(m_instMatricesInv[instId], to_float3(posAndNear));
    float3 ray_dir = mul3x3(m_instMatricesInv[instId], to_float3(dirAndFar)); // DON'T NORMALIZE IT !!!! When we transform to local space of node, ray_dir must be unnormalized!!!

    //
    //
  nanort::Ray<float> ray;

  ray.min_t = posAndNear.w;
    ray.max_t = tMax;
    ray.org[0] = ray_pos.x;
    ray.org[1] = ray_pos.y;
    ray.org[2] = ray_pos.z;
    ray.dir[0] = ray_dir.x;
    ray.dir[1] = ray_dir.y;
    ray.dir[2] = ray_dir.z;

    nanort::TriangleIntersector<float> triangle_intersector(m_verts.data() + (startSize.y * 3), m_indices.data() + startSize.x, sizeof(float) * 3);
    nanort::TriangleIntersection<float> isect;
    bool h = m_accels[geomId].Traverse(ray, triangle_intersector, &isect);
    if (h)
    {
      hit.t = isect.t;
      hit.geomId = geomId;
      hit.instId = instId;
      hit.primId = isect.prim_id;
      hit.coords[1] = isect.u;
      hit.coords[0] = isect.v;
      hit.coords[2] = 1.0f - isect.v - isect.u;
    }
  }

  return hit;
}

bool NanoRT::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  #ifdef ENABLE_METRICS
  m_stats.raysNumber++;
  #endif

  std::cout << "[NanoRT::RayQuery_AnyHit]: " 
            << "not implemeted!" << std::endl;
  return false;
}

ISceneObject* MakeNanoRT(const char* a_implName) 
  { 
  return new NanoRT();
  }

ISceneObject* MakeNanoRTExt(const char* a_implName) 
  { 
  return new NanoRTExt();
  }


///  It is here that workspaces are allocated
void NanoRTExt::CommitScene(BuildQuality a_qualityLevel)
  {
  NanoRT::CommitScene(a_qualityLevel);

  // Initialize workspace. 
  int num_threads = 1;

  ///#ifndef _DEBUG
  //num_threads = omp_get_max_threads();
  ///#endif

  m_wksp1_thr.reserve(num_threads);
  m_wksp1_thr.resize(num_threads);
  m_wksp2_thr.reserve(num_threads);
  m_wksp2_thr.resize(num_threads);

  m_w0.resize(m_accels.size());
  m_num_hits.resize(m_accels.size());
  for (size_t i = 0; i < m_accels.size(); i++)
    {
    m_w0[i].resize(((const nanortext::BVHAccelExt<float> &)m_accels[i]).GetNumNodes());
    m_num_hits[i].resize(((const nanortext::BVHAccelExt<float> &)m_accels[i]).GetNumNodes());
    for (size_t j = 0; j < m_w0[i].size(); j++)
      {
      m_w0[i][j] = 0.1;
      m_num_hits[i][j] = 0.1;
      }
    }

  // reset stats
  // 
  #ifdef ENABLE_METRICS
  nanortext::stats.Clear();
  #endif

  for (size_t i = 0; i < m_accels.size(); i++)
    m_stats.bvhTotalSize += ((const nanortext::BVHAccelExt<float> &)m_accels[i]).GetNumNodes() * sizeof(nanort::BVHNode<float>);
  m_stats.geomTotalSize = m_verts.size() * sizeof(float) + m_indices.size() * sizeof(uint32_t);
  }

CRT_Hit NanoRTExt::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
  {
  int thread_id = 0;
  ///#ifndef _DEBUG
  //thread_id = omp_get_thread_num();
  ///#endif
  nanort::BVHTraceOptions tr_opt;

  #ifdef ENABLE_METRICS
  m_stats.raysNumber++;
  ResetVarLC();
  #endif

  float3 invRayDir = SafeInverse(to_float3(dirAndFar));
  float3 rayPos = to_float3(posAndNear);

  std::vector<BoxHit>& boxMinHits = m_wksp0_thr[thread_id];
  boxMinHits.clear();
  boxMinHits.reserve(64);

  const float tMin = posAndNear.w;
  const float tMax = dirAndFar.w;

  // (1) check all instance boxes and all triangles inside box. 
  for (uint32_t i = 0; i < uint32_t(m_instBoxes.size()); i++)
    {
    const float2 minMax = RayBoxIntersection2(rayPos, invRayDir, to_float3(m_instBoxes[i].boxMin), to_float3(m_instBoxes[i].boxMax));
    if ((minMax.x <= minMax.y) && (minMax.y >= tMin) && (minMax.x <= tMax))
      boxMinHits.emplace_back(make_BoxHit(i, minMax.x));
    }

  // (2) sort all hits by hit distance to process nearest first
  //
  std::sort(boxMinHits.begin(), boxMinHits.end(), [](const BoxHit a, const BoxHit b) { return a.tHit < b.tHit; });

  // (3) process all potential hits
  //
  CRT_Hit hit;
  hit.t = tMax;
  hit.primId = uint32_t(-1);
  hit.instId = uint32_t(-1);
  hit.geomId = uint32_t(-1);

  std::vector<int>& history     = m_wksp1_thr[thread_id];
  std::vector<int>& history_aux = m_wksp2_thr[thread_id];
  history.clear();
  history_aux.clear();

  for (uint32_t boxId = 0; boxId < uint32_t(boxMinHits.size()); boxId++)
    {
    if (boxMinHits[boxId].tHit > hit.t) // already found hit that is closer than bounding box hit
      break;

    const uint32_t instId = boxMinHits[boxId].id;
    const uint32_t geomId = m_geomIdByInstId[instId];
    const uint2 startSize = m_geomOffsets[geomId];

    // transform ray with matrix to local space
    //
    float3 ray_pos = mul4x3(m_instMatricesInv[instId], to_float3(posAndNear));
    float3 ray_dir = mul3x3(m_instMatricesInv[instId], to_float3(dirAndFar)); // DON'T NORMALIZE IT !!!! When we transform to local space of node, ray_dir must be unnormalized!!!

    //
    //
    nanort::Ray<float> ray;

    ray.min_t = posAndNear.w;
    ray.max_t = tMax;
    ray.org[0] = ray_pos.x;
    ray.org[1] = ray_pos.y;
    ray.org[2] = ray_pos.z;
    ray.dir[0] = ray_dir.x;
    ray.dir[1] = ray_dir.y;
    ray.dir[2] = ray_dir.z;

    nanort::TriangleIntersector<float> triangle_intersector(m_verts.data() + (startSize.y * 3), m_indices.data() + startSize.x, sizeof(float) * 3);
    nanort::TriangleIntersection<float> isect;

    bool h = ((const nanortext::BVHAccelExt<float> &)m_accels[geomId]).Traverse(ray, triangle_intersector, &isect, tr_opt, &history_aux);
    // Increment accumuated number of hits for the geometry instance related to the ray hit
    std::vector<double>& w0 = m_w0[geomId];
    w0.resize(((const nanortext::BVHAccelExt<float>&)m_accels[geomId]).GetNumNodes());
    std::vector<real3<double>>& num_hits = m_num_hits[geomId];
    num_hits.resize(((const nanortext::BVHAccelExt<float>&)m_accels[geomId]).GetNumNodes());
    
    #ifdef ENABLE_METRICS
    int node_idx_old = 0;
    #endif

    const int nnodes = w0.size();
    for (size_t i = 0; i < history_aux.size(); i++)
      {
      const int node_idx = history_aux[i] & START_MASK;
      Assert(node_idx >= 0 && node_idx < nnodes);
      int side     = (history_aux[i] & END_MASK) >> 24;
      Assert(side >= 1 && side <= 4);
      w0[node_idx] += 1;
      if (side - 2 == -1) // the node had been queried but the ray missed its bpx
        side = 0;
      num_hits[node_idx][side - 2] += 1;
      #ifdef ENABLE_METRICS
      for (int j = 0; j < TREELET_ARR_SIZE; j++) 
      {
        if (std::abs(node_idx - node_idx_old) * sizeof(nanort::BVHNode<float>) >= size_t(treelet_sizes[j]))
          m_stats.LJC[j]++;
        const uint32_t oldCacheLineId = uint32_t(node_idx_old*sizeof(nanort::BVHNode<float>))/uint32_t(treelet_sizes[j]);
        const uint32_t newCacheLineId = uint32_t(node_idx*sizeof(nanort::BVHNode<float>))/uint32_t(treelet_sizes[j]);
        if(oldCacheLineId != newCacheLineId) 
        {
          m_stats.CMC[j]++;
          m_stats.WSS[j].insert(newCacheLineId);
        }
      }
      node_idx_old = node_idx;
      #endif
      }
    if (h)
      {
      hit.t = isect.t;
      hit.geomId = geomId;
      hit.instId = instId;
      hit.primId = isect.prim_id;
      hit.coords[1] = isect.u;
      hit.coords[0] = isect.v;
      hit.coords[2] = 1.0f - isect.v - isect.u;

      history.resize(history_aux.size());
      for (size_t j = 0; j < history_aux.size(); j++)
        history[j] = history_aux[j];
      }
    }

  return hit;
  }

/// Prints data for optimization of BVH tree, if available
int NanoRTExt::PrintTreeOptimizData(const char* fname) const
{
  if (m_w0.size() > 0)
  {
    FILE* fd = fopen(fname, "w");
    for (size_t ig = 0; ig < m_accels.size(); ig++)
      {
      fprintf(fd, "Geometry ID = %10i\n", int(ig));
      for (size_t i = 0; i < m_w0[ig].size(); i++)
        fprintf(fd, "%20i %15.3f %15.6f %15.3f\n", int(i), m_w0[ig][i],
          ((const nanortext::BVHAccelExt<float> &)m_accels[ig]).NodeArea(i),
          m_w0[ig][i] / ((const nanortext::BVHAccelExt<float> &)m_accels[ig]).NodeArea(i));
      fprintf(fd, "\n\n\n\n\n");
      }
    fclose(fd);
  }
  return 0;
}

/// Optimize tree, if has data for it
int NanoRTExt::OptimizeBVH()
  {
  for (size_t ig = 0; ig < m_accels.size(); ig++)
    {
    std::vector<double>& w0 = m_w0[ig];
    std::vector< real3<double> >& num_hits = m_num_hits[ig];
    if (w0.size() == 0)
      continue; // no data

#if 0
    if (((nanortext::BVHAccelExt<float> &)m_accels[ig]).ImproveSmallNodesProb(0, num_hits) != SUCCESS)
      {
      Assert(false);
      return FAILURE;
      }
#endif

    //std::vector<double> w0_new;
    //w0_new.resize(w0.size());
    for (size_t i = 0; i < w0.size(); i++)
      w0[i] = num_hits[i][0] + num_hits[i][1] + num_hits[i][2];

    PrintTreeOptimizData("w0_new.dat");

    if (((nanortext::BVHAccelExt<float> &)m_accels[ig]).ReArrangeByClusters(w0, 0.9) != SUCCESS)
      {
      Assert(false);
      return FAILURE;
      }
    }
  return SUCCESS;
  }

/// Get statistics
MetricStats NanoRTExt::GetStats()
{
  #ifdef ENABLE_METRICS
  m_stats.BLB = nanortext::stats.BLB;
  m_stats.LC  = nanortext::stats.LC;
  m_stats.LC2 = nanortext::stats.LC2;
  for (int i = 0; i < TREELET_ARR_SIZE; i++)
    m_stats.LJC[i] = nanortext::stats.LJC[i];
  m_stats.NC  = nanortext::stats.NC;
  m_stats.SBL = nanortext::stats.SBL;
  m_stats.SOC = nanortext::stats.SOC;
  m_stats.TC  = nanortext::stats.TC;
  #endif
  return NanoRT::GetStats();
}