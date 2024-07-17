#include <iostream>
#include <queue>
#include <stack>
#include <limits>
#include <tuple>
#include <cmath>

#include "cbvh_core.h"
#include "cbvh_internal.h"

using LiteMath::float3;
using LiteMath::float4;
using LiteMath::float4x4;
using LiteMath::uint;

#define M_PI		3.14159265358979323846
#define TASKING_INTERNAL
#include <kernels/bvh/bvh.h>
#include <kernels/geometry/trianglev.h>
#include <embree3/rtcore.h>
#include <embree3/rtcore_ray.h>


static void customBoundsFunc(const struct RTCBoundsFunctionArguments* args)
{
  const cbvh::BVHNode* customObjects = (const cbvh::BVHNode*) args->geometryUserPtr;
  RTCBounds* bounds_o = args->bounds_o;
  const cbvh::BVHNode& obj = customObjects[args->primID];
  bounds_o->lower_x = obj.boxMin.x;
  bounds_o->lower_y = obj.boxMin.y;
  bounds_o->lower_z = obj.boxMin.z;
  bounds_o->upper_x = obj.boxMax.x;
  bounds_o->upper_y = obj.boxMax.y;
  bounds_o->upper_z = obj.boxMax.z;
}

static void customIntersectFunc(const RTCIntersectFunctionNArguments* args)
{
  // for testing purposes
}

static void customOccludedFunc(const RTCOccludedFunctionNArguments* args)
{
  // for testing purposes
}

inline void createCustomGeometryObject (RTCScene scene, RTCDevice device, const cbvh::BVHNode& in_bvhNode)
{
  RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
  //in_bvhNode.escapeIndex = rtcAttachGeometry(scene, geom);
  rtcAttachGeometry(scene, geom);  

  rtcSetGeometryUserPrimitiveCount(geom, 1);
  rtcSetGeometryUserData(geom, (void*)&in_bvhNode);
  rtcSetGeometryBoundsFunction(geom, customBoundsFunc, nullptr);

  rtcSetGeometryIntersectFunction(geom, customIntersectFunc);
  rtcSetGeometryOccludedFunction(geom, customOccludedFunc);

  rtcCommitGeometry(geom);

  rtcReleaseGeometry(geom);
  return;
}


namespace cbvh
{
  void GetBVH4(RTCScene& in_rtcScene, const uint* a_indices, size_t a_indexNum, cbvh::BVHTree& res);
};

namespace embree
{
  /* error reporting function */
  void error_handler(void* userPtr, const RTCError code, const char* str)
  {
    if (code == RTC_ERROR_NONE) 
      return;
    
    printf("Embree: ");
    switch (code) {
    case RTC_ERROR_UNKNOWN          : printf("RTC_ERROR_UNKNOWN"); break;
    case RTC_ERROR_INVALID_ARGUMENT : printf("RTC_ERROR_INVALID_ARGUMENT"); break;
    case RTC_ERROR_INVALID_OPERATION: printf("RTC_ERROR_INVALID_OPERATION"); break;
    case RTC_ERROR_OUT_OF_MEMORY    : printf("RTC_ERROR_OUT_OF_MEMORY"); break;
    case RTC_ERROR_UNSUPPORTED_CPU  : printf("RTC_ERROR_UNSUPPORTED_CPU"); break;
    case RTC_ERROR_CANCELLED        : printf("RTC_ERROR_CANCELLED"); break;
    default                         : printf("invalid error code"); break;
    }
    if (str) { 
      printf(" ("); 
      while (*str) putchar(*str++); 
      printf(")\n"); 
    }
    exit(1);
  }


};

#ifdef WIN32
namespace embree // forgotten piece of embree shit to make this works on windows
{
  TrueTy True;
  FalseTy False;
  ZeroTy zero;
  OneTy one;
  NegInfTy neg_inf;
  PosInfTy inf;
  PosInfTy pos_inf;
  NaNTy nan;
  UlpTy ulp;
  PiTy pi;
  OneOverPiTy one_over_pi;
  TwoPiTy two_pi;
  OneOverTwoPiTy one_over_two_pi;
  FourPiTy four_pi;
  OneOverFourPiTy one_over_four_pi;
  StepTy step;
  ReverseStepTy reverse_step;
  EmptyTy empty;
  UndefinedTy undefined;
}
#endif

RTCDevice g_intelPieceOfShit = nullptr;
RTCScene  g_scene = nullptr;

/**
  \brief Initializes Embree scene and device instances
*/
void cbvh::Init()
{
  g_intelPieceOfShit = rtcNewDevice("threads=1, tri_accel=bvh4.triangle4v");
  g_scene            = rtcNewScene(g_intelPieceOfShit);
}

/**
  \brief Releases Embree scene and device instances
*/
void cbvh::Destroy()
{
  rtcReleaseScene (g_scene);            g_scene = nullptr;
  rtcReleaseDevice(g_intelPieceOfShit); g_intelPieceOfShit = nullptr;
}

/**
  \brief
*/
struct NodePair
{
  NodePair() : theirs(embree::BVH4::NodeRef()), oursOffs(0xFFFFFFFF) { }
  NodePair(embree::BVH4::NodeRef a, uint b) : theirs(a), oursOffs(b) {}
  embree::BVH4::NodeRef theirs;
  uint oursOffs;
};

/*
 * Bottom level Tree Depth-First Traversal
 */
struct DFSParamsBottomLevel {
    DFSParamsBottomLevel(const embree::BVH4 *a_bvhTheirs,
                         const std::vector<cbvh::Interval> &a_primsBegEnd, const std::vector<uint> &a_primIdBFS,
                         const uint *a_indices, size_t a_indexNum,
                         cbvh::BVHTree &a_bvhOurs) : m_bvhTheirs(a_bvhTheirs),
                                          m_primsBegEnd(a_primsBegEnd), m_primIdBFS(a_primIdBFS),
                                          m_indices(a_indices), m_indexNum(a_indexNum), m_bvhOurs(a_bvhOurs) {}

    size_t BFSTraversalRec(NodePair currNode) {
      const cbvh::BVHNode &nodeOurs = m_bvhOurs.nodes[currNode.oursOffs];
      // if the current node is not a leaf
      if (currNode.theirs.isAlignedNode()) {
        assert(!cbvh::IsLeaf(nodeOurs));

        for (int i = 0; i < 4; i++) 
        {
          const size_t start = m_bvhOurs.indicesReordered.size()/3;
          const size_t end   = BFSTraversalRec(NodePair(currNode.theirs.alignedNode()->child(i), nodeOurs.leftOffset + i));

          // write down intervals, basically taking the sum from the bottom levels of tree
          m_bvhOurs.intervals[nodeOurs.leftOffset + i].start = uint(start);
          m_bvhOurs.intervals[nodeOurs.leftOffset + i].count = uint(end - start);
        }
      } else {
        assert(cbvh::IsLeaf(nodeOurs) || cbvh::IsEmpty(nodeOurs));

        const cbvh::Interval begEnd = m_primsBegEnd[currNode.oursOffs];

        for (int tri = begEnd.start; tri < (begEnd.start + begEnd.count); tri++) {
          const uint primId = m_primIdBFS[tri];

          // write down reordered indices
          m_bvhOurs.indicesReordered.push_back(m_indices[primId * 3 + 0]);
          m_bvhOurs.indicesReordered.push_back(m_indices[primId * 3 + 1]);
          m_bvhOurs.indicesReordered.push_back(m_indices[primId * 3 + 2]);
        }

      } // end if

      return m_bvhOurs.indicesReordered.size()/3;
    }

    const embree::BVH4 *m_bvhTheirs;
    const std::vector<cbvh::Interval> &m_primsBegEnd;
    const std::vector<uint> &m_primIdBFS;
    const uint *m_indices;
    size_t m_indexNum;
    cbvh::BVHTree &m_bvhOurs;
};

void SimultaneousDFSTraversalBottomLevel(const embree::BVH4 *a_bvhTheirs,
                                         const std::vector<cbvh::Interval> &a_primsBegEnd, const std::vector<uint> &a_primIdBFS,
                                         const uint *a_indices, size_t a_indexNum,
                                         cbvh::BVHTree &a_bvhOurs) {

  a_bvhOurs.indicesReordered.resize(0);
  a_bvhOurs.intervals.resize(a_bvhOurs.nodes.size());

  DFSParamsBottomLevel params(a_bvhTheirs, a_primsBegEnd, a_primIdBFS, a_indices, a_indexNum, a_bvhOurs);

  params.BFSTraversalRec(NodePair(a_bvhTheirs->root, 0));

  a_bvhOurs.intervals[0] = cbvh::Interval(0, a_indexNum / 3);
}

/*
 * Top level Tree Depth-First Traversal
 */
struct DFSParamsTopLevel {
    DFSParamsTopLevel(const embree::BVH4 *a_bvhTheirs,
                      const std::vector<cbvh::Interval> &a_primsBegEnd, const std::vector<uint> &a_primIdBFS,
                      const uint *a_indices, size_t a_indexNum,
                      cbvh::BVHTree &a_bvhOurs) : m_bvhTheirs(a_bvhTheirs),
                                                  m_primsBegEnd(a_primsBegEnd), m_primIdBFS(a_primIdBFS),
                                                  m_indices(a_indices), m_indexNum(a_indexNum), m_bvhOurs(a_bvhOurs) {}

    size_t BFSTraversalRec(NodePair currNode) {
      const cbvh::BVHNode &nodeOurs = m_bvhOurs.nodes[currNode.oursOffs];
      // if the current node is not a leaf
      if (currNode.theirs.isAlignedNode()) {
        assert(!cbvh::IsLeaf(nodeOurs));

        for (int i = 0; i < 4; i++)
        {
          const size_t start = m_bvhOurs.indicesReordered.size();
          const size_t end   = BFSTraversalRec(NodePair(currNode.theirs.alignedNode()->child(i), nodeOurs.leftOffset + i));

          // write down intervals, basically taking the sum from the bottom levels of tree
          m_bvhOurs.intervals[nodeOurs.leftOffset + i].start = uint(start);
          m_bvhOurs.intervals[nodeOurs.leftOffset + i].count = uint(end - start);
        }
      } else {
        //assert(cbvh::IsLeaf(nodeOurs));
        if (!cbvh::IsEmpty(nodeOurs)) {
          const cbvh::Interval begEnd = m_primsBegEnd[currNode.oursOffs];
          const uint primId = m_primIdBFS[begEnd.start];

          // write down reordered indices
          m_bvhOurs.indicesReordered.push_back(m_indices[primId]);
        }
      } // end if

      return m_bvhOurs.indicesReordered.size();
    }

    const embree::BVH4 *m_bvhTheirs;
    const std::vector<cbvh::Interval> &m_primsBegEnd;
    const std::vector<uint> &m_primIdBFS;
    const uint *m_indices;
    size_t m_indexNum;
    cbvh::BVHTree &m_bvhOurs;
};

void SimultaneousDFSTraversalTopLevel(const embree::BVH4 *a_bvhTheirs,
                              const std::vector<cbvh::Interval> &a_primsBegEnd, const std::vector<uint> &a_primIdBFS,
                              const uint *a_indices, size_t a_indexNum,
                              cbvh::BVHTree &a_bvhOurs) {

  a_bvhOurs.indicesReordered.resize(0);
  a_bvhOurs.intervals.resize(a_bvhOurs.nodes.size());

  DFSParamsTopLevel params(a_bvhTheirs, a_primsBegEnd, a_primIdBFS, a_indices, a_indexNum, a_bvhOurs);

  params.BFSTraversalRec(NodePair(a_bvhTheirs->root, 0));

  a_bvhOurs.intervals[0] = cbvh::Interval(0, a_indexNum);
}

// A temporary structure, used in GetBVH4 to pass the node with current depth and escape index
struct TempNode{
    TempNode(embree::BVH4::NodeRef in_node, int in_escapeIndex, int in_currentDepth)
      : node(in_node), escapeIndex(in_escapeIndex), currentDepth(in_currentDepth){}
    embree::BVH4::NodeRef node;
    int escapeIndex;
    int currentDepth;
};

void cbvh::GetBVH4(RTCScene& in_rtcScene, const uint* a_indices, size_t a_indexNum, cbvh::BVHTree& res)
{
  embree::BVH4 *bvh4 = nullptr;

  // Should be used for bottom level bvh, so TY_BVH4 it is
  embree::AccelData *accel = ((embree::Accel *) in_rtcScene)->intersectors.ptr;
  if (accel->type == embree::AccelData::TY_BVH4)
    bvh4 = (embree::BVH4 *) accel;
  else
    return; // #TODO: need error code

  // Queue with tuples of the form: { Embree BVH node; Escape index; Current depth }
  std::queue<TempNode> treeQueue;

  const size_t nodesApproxNum = a_indexNum/2;
  res.nodes.reserve(nodesApproxNum);
  res.intervals.reserve(nodesApproxNum);
  res.indicesReordered.reserve(a_indexNum);

  std::vector<uint> primIdBFS;
  std::vector<cbvh::Interval> primsBegEnd;

  primIdBFS.reserve(a_indexNum/3 + 10);
  primsBegEnd.reserve(nodesApproxNum);

  //cvex4::set_ftz();

  int currentLevelDepth = 0;
  int currentLevelStartIndex = 0;

  // In the first pass we traverse the tree breath-first using the queue
  treeQueue.push(TempNode(bvh4->root, uint(-2), 0));
  while (treeQueue.size() != 0) 
  {
    // check if the current node has jumped to the next depth level
    if(treeQueue.front().currentDepth > currentLevelDepth)
    {
      Interval currentLevelDepthInterval(currentLevelStartIndex, res.nodes.size() - currentLevelStartIndex);
      // push the interval of the previous depth level
      res.depthRanges.push_back(currentLevelDepthInterval);
      currentLevelDepth += 1;
      currentLevelStartIndex = res.nodes.size();
    }
    // if the current node is really a node and not a leaf
    if(treeQueue.front().node.isAlignedNode())
    { 
      const int escapeIndex = treeQueue.front().escapeIndex;
      const int bufferSize  = res.nodes.size() + treeQueue.size();
      
      const auto* node      = treeQueue.front().node.alignedNode();
      const auto& boxMin    = node->bounds().lower;
      const auto& boxMax    = node->bounds().upper;

      // push the current node to the resulting bvhtree
      res.nodes.push_back(BVHNode(float3(boxMin.x, boxMin.y, boxMin.z), 
                                  float3(boxMax.x, boxMax.y, boxMax.z), 
                                  bufferSize, escapeIndex));
      primsBegEnd.push_back(cbvh::Interval());

      // push the child nodes into the queue
      for (size_t i = 0; i < 3; i++) 
        treeQueue.push(TempNode(node->child(i), bufferSize + i + 1, treeQueue.front().currentDepth + 1));
      
      treeQueue.push(TempNode(node->child(3), escapeIndex, treeQueue.front().currentDepth + 1));
      treeQueue.pop();
    } 
    else // if the current node is the leaf
    {
      size_t num;
      const embree::Triangle4v *tri = (const embree::Triangle4v *)treeQueue.front().node.leaf(num);

      float4 bMin = float4( +std::numeric_limits<float>::infinity() );
      float4 bMax = float4( -std::numeric_limits<float>::infinity() ); 
      
      size_t primStart = primIdBFS.size();
      for (size_t i = 0; i < num; i++) // #TODO: opt min/max code for xxxx yyyy zzzz embree data layout
      {
        for (size_t j = 0; j < tri[i].size(); j++) 
        {
          const auto primId = tri[i].primID(j);
          primIdBFS.push_back(primId);

          float4 A = float4{ tri[i].v0.x[j], tri[i].v0.y[j], tri[i].v0.z[j], 0.0f }; // not optimal but working approach
          float4 B = float4{ tri[i].v1.x[j], tri[i].v1.y[j], tri[i].v1.z[j], 0.0f }; // not optimal but working approach
          float4 C = float4{ tri[i].v2.x[j], tri[i].v2.y[j], tri[i].v2.z[j], 0.0f }; // not optimal but working approach

          bMin = min(min(bMin,A), min(B,C));
          bMax = max(max(bMax,A), max(B,C));
        } 
      }
      size_t primEnd = primIdBFS.size();

      const float4 boxMin   = bMin;
      const float4 boxMax   = bMax;
      const int escapeIndex = treeQueue.front().escapeIndex;

      unsigned int offset_flag = 0;
      if (num > 0 )
        offset_flag = 0xFFFFFFFF;
      else{
        offset_flag = 0xFFFFFFFD;
      }

      res.nodes.push_back(BVHNode(float3(boxMin.x, boxMin.y, boxMin.z),
                                  float3(boxMax.x, boxMax.y, boxMax.z),
                                  offset_flag, escapeIndex));
      
      primsBegEnd.push_back( cbvh::Interval(uint(primStart), uint(primEnd - primStart)) );
      treeQueue.pop();
    }
  }
  // Add the very last level of depth
  Interval currentLevelDepthInterval(currentLevelStartIndex, res.nodes.size() - currentLevelStartIndex);
  // push the interval of the previous depth level
  res.depthRanges.push_back(currentLevelDepthInterval);

  // Second pass with depth-first traversal to form triangle indices
  SimultaneousDFSTraversalBottomLevel(bvh4, primsBegEnd, primIdBFS, a_indices, a_indexNum, res);

  return;
}

cbvh::BVHTree cbvh::BuildBVH4(const float4* a_vertices, size_t a_vertNum, const uint* a_indices, size_t a_indexNum) {
  assert(g_scene != nullptr);
  assert(g_intelPieceOfShit != nullptr);

  /// Loading mesh to Embree
  RTCGeometry mesh = rtcNewGeometry(g_intelPieceOfShit, RTC_GEOMETRY_TYPE_TRIANGLE);

  float4 *vertices = (float4 *) rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
                                                        sizeof(float4), a_vertNum);
  uint *triangles = (uint *) rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(uint) * 3,
                                                     a_indexNum / 3);

  memcpy(vertices, a_vertices, sizeof(float4) * a_vertNum);
  memcpy(triangles, a_indices, sizeof(uint) * a_indexNum);

  rtcCommitGeometry(mesh);

  // Resulting bvh tree
  cbvh::BVHTree res;
  res.geomID = rtcAttachGeometry(g_scene, mesh);
  rtcReleaseGeometry(mesh);

  rtcSetSceneBuildQuality(g_scene, RTC_BUILD_QUALITY_HIGH);
  rtcCommitScene(g_scene);

  GetBVH4(g_scene, a_indices, a_indexNum,res);

  return res;
}

cbvh::BVHTree cbvh::BuildBVH4(const BVHNode* a_nodes, size_t a_objNum)
{
  std::vector<uint> objIndices(a_objNum);
  for(size_t i=0;i<a_objNum;i++)
    objIndices[i] = a_nodes[i].leftOffset;

  const uint* a_indices = objIndices.data(); 

  // Create and commit to Embree custom geometry
  for(int i = 0; i < a_objNum; ++i) {
    createCustomGeometryObject(g_scene, g_intelPieceOfShit, a_nodes[i]);
  }

  rtcCommitScene(g_scene);

  //PrintBVH(g_scene);

  cbvh::BVHTree res;

  embree::BVH4 *bvh4 = nullptr;

  // TY_BVH4 is always used here since we don't load meshes in the same tree but represent them with bbox
  embree::AccelData *accel = ((embree::Accel *) g_scene)->intersectors.ptr;
  if (accel->type == embree::AccelData::TY_BVH4)
    bvh4 = (embree::BVH4 *) accel;
  else
    return res; // #TODO: need error code

  // Queue with tuples of the form: { Embree BVH node; Escape index; Current depth }
  std::queue<TempNode> treeQueue;

  const size_t nodesApproxNum = a_objNum * 1.5;
  res.nodes.reserve(nodesApproxNum);
  res.intervals.reserve(nodesApproxNum);
  res.indicesReordered.reserve(a_objNum);


  std::vector<uint> primIdBFS;
  std::vector<cbvh::Interval> primsBegEnd;

  primIdBFS.reserve(a_objNum + 10);
  primsBegEnd.reserve(nodesApproxNum);

  //cvex4::set_ftz();

  int currentLevelDepth = 0;
  int currentLevelStartIndex = 0;

  // In the first pass we traverse the tree breath-first using the queue
  treeQueue.push(TempNode(bvh4->root, uint(-2), 0));
  while (treeQueue.size() != 0)
  {
    // check if the current node has jumped to the next depth level
    if(treeQueue.front().currentDepth > currentLevelDepth)
    {
      Interval currentLevelDepthInterval(currentLevelStartIndex, res.nodes.size() - currentLevelStartIndex);
      // push the interval of the previous depth level
      res.depthRanges.push_back(currentLevelDepthInterval);
      currentLevelDepth += 1;
      currentLevelStartIndex = res.nodes.size();
    }
    // if the current node is really a node and not a leaf
    if(treeQueue.front().node.isAlignedNode())
    {
      const int escapeIndex = treeQueue.front().escapeIndex;
      const int bufferSize  = res.nodes.size() + treeQueue.size();

      const auto* node      = treeQueue.front().node.alignedNode();
      const auto& boxMin    = node->bounds().lower;
      const auto& boxMax    = node->bounds().upper;

      // push the current node to the resulting bvhtree
      res.nodes.push_back(BVHNode(float3(boxMin.x, boxMin.y, boxMin.z),
                                  float3(boxMax.x, boxMax.y, boxMax.z),
                                  bufferSize, escapeIndex));

      primsBegEnd.push_back(cbvh::Interval());

      // push the child nodes into the queue
      for (size_t i = 0; i < 3; i++)
        treeQueue.push(TempNode(node->child(i), bufferSize + i + 1, treeQueue.front().currentDepth + 1));

      treeQueue.push(TempNode(node->child(3), escapeIndex, treeQueue.front().currentDepth + 1));
      treeQueue.pop();
    }
    else // if the current node is the leaf
    {
      size_t num;
      const int* nodeType = (const int*)treeQueue.front().node.leaf(num);
      //const embree::Triangle4v *tri = (const embree::Triangle4v *)treeQueue.front().node.leaf(num);
      auto pinf = +std::numeric_limits<float>::infinity();
      auto ninf = -std::numeric_limits<float>::infinity();

      size_t primStart = primIdBFS.size();

      // In case the node is empty
      if(nodeType == NULL) {

        primIdBFS.push_back(0xFFFFFFFD);
        res.nodes.push_back(BVHNode(float3(pinf, pinf, pinf),
                                    float3(ninf, ninf, ninf),
                                    0xFFFFFFFD, treeQueue.front().escapeIndex));
        primsBegEnd.push_back(cbvh::Interval(uint(primStart), 0));

        treeQueue.pop();
      }
      else {
        primIdBFS.push_back(*nodeType);

        const float3 boxMin = a_nodes[*nodeType].boxMin;
        const float3 boxMax = a_nodes[*nodeType].boxMax;

        res.nodes.push_back(BVHNode(float3(boxMin.x, boxMin.y, boxMin.z),
                                    float3(boxMax.x, boxMax.y, boxMax.z),
                                    0xFFFFFFFF, treeQueue.front().escapeIndex));

        primsBegEnd.push_back(cbvh::Interval(uint(primStart), 1));

        treeQueue.pop();
      }
    }
  }
  // Add the very last level of depth
  Interval currentLevelDepthInterval(currentLevelStartIndex, res.nodes.size() - currentLevelStartIndex);
  // push the interval of the previous depth level
  res.depthRanges.push_back(currentLevelDepthInterval);

  SimultaneousDFSTraversalTopLevel(bvh4, primsBegEnd, primIdBFS, a_indices, a_objNum,res);

  return res;
}
