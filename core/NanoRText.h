#ifndef NANORTEXT_H_
#define NANORTEXT_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <vector>

#include "raytrace_common.h"

static const int SUCCESS = 0;
static const int FAILURE = -1;

void Assert(bool expr)
  {
  assert(expr);
  }

#include "nanort/nanort.h"

using namespace nanort;

namespace nanortext {

#ifdef ENABLE_METRICS
  struct Metrics
    {
    uint32_t NC;  ///< Nodes Count
    uint32_t LC;  ///< Leaves Count
    uint32_t LC2;  ///< Leaves Count2
    uint32_t TC;  ///< Triangles Count
    uint32_t LJC[TREELET_ARR_SIZE]; ///< Long Jumps Count -- number of jumps far then treelet size during traversal;
                  ///< please note that we also account jumps on stack.pop in this metric, howerev we can discuss this.
    uint32_t BLB; ///< Bus Load in bytes
    uint32_t SOC; ///< Stack Operations Count (both push and pop)
    uint32_t SBL; ///< Stack Bytes Load (both push and pop)

    void Clear() { Metrics zerows = {}; *this = zerows; }
    };

static Metrics stats;
#endif

/// @brief Extended Bounding Volume Hierarchy acceleration which provides tree 
///                 re-arrangement in memory by clusterization of nodes using
///                 "weight" of nodes as criterion.
///
/// @tparam T real value type(float or double).
///
template <typename T>
class BVHAccelExt : public nanort::BVHAccel<T> {
 public:

  /// @brief Traverse into BVH along ray and find closest hit point & primitive if
  template <class I, class H>
  bool Traverse(const Ray<T> &ray, const I &intersector, H *isect,
                const BVHTraceOptions &options = BVHTraceOptions(),
                std::vector<int>* history = NULL) const;

  /// The number of nodes
  int GetNumNodes() const 
    {
    return BVHAccel<T>::nodes_.size();
    };

  /// The number of nodes in the subtree
  int NumSubTreeNodes(int root_node_idx) const;

  /// Area of the given (i-th) node. Sun of area of 3 walls of the bounding box
  double NodeArea(int i) const
    {
    const BVHNode<T>& node = BVHAccel<T>::nodes_[i];
    double area;
    area = (node.bmax[0] - node.bmin[0]) * (node.bmax[1] - node.bmin[1])
         + (node.bmax[1] - node.bmin[1]) * (node.bmax[2] - node.bmin[2])
         + (node.bmax[0] - node.bmin[0]) * (node.bmax[2] - node.bmin[2]);
    return area;
    };

  protected:
  /// Determines which nodes (of a sub-tree) constitute its root cluster
  /// and which (those "below" the cluster in the tree hierarchy) will be
  /// the roots of child clusters.
  ///
  /// Both arrays of nodes are sorted so that the "criterion" w0 (can be arbitrary
  /// but aimed mainly for probability of ray hitting the node box) is descending.
  ///
  /// @param[in] root_node_idx - index of the root of this subtree
  /// @param[in] w0            - weights of nodes ("criterion")
  /// @param[in] threshold     - threshold to select nodes for the root cluster
  /// @param[in] max_num_nodes - max allowed number of cluster nodes (but it can contain less)
  /// @param[out] rcluster_idx - indices of the nodes that will form the root cluster
  ///                            Sorted according to criterion
  /// @param[out] child_cluster_root_idx - indices of the nodes that will become roots of
  ///                                      the child clusters.
  ///                                      Sorted according to criterion
  /// @return SUCCESS/FAILURE
  int RootCluster4Subtree(int root_node_idx,
                          const std::vector<double> &w0,
                          double threshold,
                          int max_num_nodes,
                          std::vector<int> &rcluster_idx,
                          std::vector<int> &child_cluster_root_idx) const;

  /// Suggests ordering of the nodes of a subtree "by clusters". First comes the
  /// "root cluster" then "child clusters". Each of the latter, in turn, consists
  /// of its root cluster, past which come child clusters, of which... etc.
  ///
  /// If a cluster is NOT further subdivided, its nodes are ordered so that their
  /// criterion "w0" is descending.
  /// If it is further clusterized, then first comes its root cluster, then child
  /// clusters in such an order that first comes the one whose root has the largest
  /// "w0", then the one whose root has the next "w0" and so on.
  ///
  /// The order of nodes inside each sub-cluster (root or child) obeys the same
  /// obey the same rule, i.e. if this sub-cluster is further subdivided, then... else
  /// ordered so that their "w0" is descending.
  ///
  /// @param[in] root_node_idx - index of the root of this subtree
  /// @param[in] w0            - weights of nodes ("criterion")
  /// @param[in] threshold     - threshold to select nodes for the root cluster
  /// @param[in] max_num_nodes - max allowed number of cluster nodes (but it can contain less)
  /// @param[out] indices      - indices of ALL nodes of the sub-tree in the new order
  ///
  /// @return SUCCESS/FAILURE
  int ClusterizeSubtree(int root_node_idx,
                        const std::vector<double> &w0,
                        double threshold,
                        std::vector<int> &indices) const;

  public:
  /// Changes the order of the tree nodes in memory WHLE RETAINING THE TREE TOPOLOGY
  /// and updates the links to the children kept in the nodes accordingly.
  ///
  /// New order of nodes is that they are stored "cluster by cluster", each cluster,
  /// in turn, is usually further subdivided by sub-clusters etc. But if it is too small
  /// it is not further subdivided.
  ///
  /// If a cluster is NOT further subdivided, its nodes are ordered so that their
  /// criterion "w0" is descending.
  /// If it is further clusterized, then first comes its root cluster, then child
  /// clusters in such an order that first comes the one whose root has the largest
  /// "w0", then the one whose root has the next "w0" and so on.
  ///
  /// The order of nodes inside each sub-cluster (root or child) obeys the same
  /// obey the same rule, i.e. if this sub-cluster is further subdivided, then... else
  /// ordered so that their "w0" is descending.
  ///
  /// @param[in] w0            - weights of nodes ("criterion")
  /// @param[in] threshold     - threshold to select nodes for the root cluster
  /// @return SUCCESS/FAILURE
  int ReArrangeByClusters(const std::vector<double> &w0, double threshold);

  /// Traverse subtree and for small-size nodes calculate frequency of hits by scaling that of parent
  int ImproveSmallNodesProb(int root_node_idx, std::vector< real3<double> > &num_hits) const;
};


/// Element of array to sort in a DESCENDING order. Keeps old index.
struct FloatAndInt
  {
  float x;
  int i;

  //  FloatAndInt() {};
  FloatAndInt(float y, int j)
    {
    x = y;
    i = j;
    };

  /// ... and its comparator for sort()
  static bool greater(const FloatAndInt& a, const FloatAndInt& b)
    {
    return a.x > b.x;
    }
  };

/// @return 0 in case ray missed the box, otherwise 1, 2, 3 if the wall hit is orthogonal to X,Y,Z
inline int IntersectRayAABBExt(float* tminOut,  // [out]
                               float* tmaxOut,  // [out]
                               float min_t, float max_t,
                               const float bmin[3], const float bmax[3],
                               const real3<float>& ray_org,
                               const real3<float>& ray_inv_dir,
                               int ray_dir_sign[3])
  {
  float tmin, tmax;

  const float min_x = ray_dir_sign[0] ? bmax[0] : bmin[0];
  const float min_y = ray_dir_sign[1] ? bmax[1] : bmin[1];
  const float min_z = ray_dir_sign[2] ? bmax[2] : bmin[2];
  const float max_x = ray_dir_sign[0] ? bmin[0] : bmax[0];
  const float max_y = ray_dir_sign[1] ? bmin[1] : bmax[1];
  const float max_z = ray_dir_sign[2] ? bmin[2] : bmax[2];

  // X
  const float tmin_x = (min_x - ray_org[0]) * ray_inv_dir[0];
  // MaxMult robust BVH traversal(up to 4 ulp).
  // 1.0000000000000004 for double precision.
  const float tmax_x = (max_x - ray_org[0]) * ray_inv_dir[0] * 1.00000024f;

  // Y
  const float tmin_y = (min_y - ray_org[1]) * ray_inv_dir[1];
  const float tmax_y = (max_y - ray_org[1]) * ray_inv_dir[1] * 1.00000024f;

  // Z
  const float tmin_z = (min_z - ray_org[2]) * ray_inv_dir[2];
  const float tmax_z = (max_z - ray_org[2]) * ray_inv_dir[2] * 1.00000024f;

  tmax = safemin(tmax_z, safemin(tmax_y, safemin(tmax_x, max_t)));

  //  tmin = safemax(tmin_z, safemax(tmin_y, safemax(tmin_x, min_t)));
  int idx = 0;
  if (tmin_x > tmin_y)
    {
    if (tmin_z > tmin_x)
      {
      idx = 3;
      tmin = safemax(tmin_z, min_t);
      }
    else
      {
      idx = 1;
      tmin = safemax(tmin_x, min_t);
      }
    }
  else
    {
    if (tmin_z > tmin_y)
      {
      idx = 3;
      tmin = safemax(tmin_z, min_t);
      }
    else
      {
      idx = 2;
      tmin = safemax(tmin_y, min_t);
      }
    }


  if (tmin <= tmax)
    {
    (*tminOut) = tmin;
    (*tmaxOut) = tmax;
    return idx;
    }

  return 0;  // no hit
  }

/// @return 0 in case ray missed the box, otherwise 1, 2, 3 if the wall hit is orthogonal to X,Y,Z
inline int IntersectRayAABBExt(double* tminOut,  // [out]
                               double* tmaxOut,  // [out]
                               double min_t, double max_t,
                               const double bmin[3], const double bmax[3],
                               const real3<double> &ray_org,
                               const real3<double> &ray_inv_dir,
                               int ray_dir_sign[3]) 
  {
  double tmin, tmax;

  const double min_x = ray_dir_sign[0] ? bmax[0] : bmin[0];
  const double min_y = ray_dir_sign[1] ? bmax[1] : bmin[1];
  const double min_z = ray_dir_sign[2] ? bmax[2] : bmin[2];
  const double max_x = ray_dir_sign[0] ? bmin[0] : bmax[0];
  const double max_y = ray_dir_sign[1] ? bmin[1] : bmax[1];
  const double max_z = ray_dir_sign[2] ? bmin[2] : bmax[2];

  // X
  const double tmin_x = (min_x - ray_org[0]) * ray_inv_dir[0];
  // MaxMult robust BVH traversal(up to 4 ulp).
  const double tmax_x = (max_x - ray_org[0]) * ray_inv_dir[0] * 1.0000000000000004;

  // Y
  const double tmin_y = (min_y - ray_org[1]) * ray_inv_dir[1];
  const double tmax_y = (max_y - ray_org[1]) * ray_inv_dir[1] * 1.0000000000000004;

  // Z
  const double tmin_z = (min_z - ray_org[2]) * ray_inv_dir[2];
  const double tmax_z = (max_z - ray_org[2]) * ray_inv_dir[2] * 1.0000000000000004;

  tmax = safemin(tmax_z, safemin(tmax_y, safemin(tmax_x, max_t)));

  //  tmin = safemax(tmin_z, safemax(tmin_y, safemax(tmin_x, min_t)));
  int idx = 0;
  if (tmin_x > tmin_y)
    {
    if (tmin_z > tmin_x)
      {
      idx = 3;
      tmin = safemax(tmin_z, min_t);
      }
    else
      {
      idx = 1;
      tmin = safemax(tmin_x, min_t);
      }
    }
  else
    {
    if (tmin_z > tmin_y)
      {
      idx = 3;
      tmin = safemax(tmin_z, min_t);
      }
    else
      {
      idx = 2;
      tmin = safemax(tmin_y, min_t);
      }
    }


  if (tmin <= tmax) 
    {
    (*tminOut) = tmin;
    (*tmaxOut) = tmax;
    return idx;
    }

  return 0;  // no hit
  }



/// @brief Traverse into BVH along ray and find closest hit point & primitive if
/// found
///
/// @tparam I Intersector class
/// @tparam H Hit class
///
/// @param[in] ray Input ray
/// @param[in] intersector Intersector object. This object is called for each possible intersection of ray and BVH during traversal.
/// @param[out] isect Intersection point information(filled when closest hit point was found)
/// @param[in] options Traversal options.
/// @param[out] history - IGNORED, if NULL. Otherwise keeps indices of the
///                       nodes tested, in the order of traversal.
///
/// @return true if the closest hit point found.
///
template <typename T>
template <class I, class H>
bool BVHAccelExt<T>::Traverse(const Ray<T>& ray, const I& intersector, H* isect,
  const BVHTraceOptions& options, std::vector<int>* history) const {
  const int kMaxStackDepth = 512;
  (void)kMaxStackDepth;

  T hit_t = ray.max_t;

  int node_stack_index = 0;
  unsigned int node_stack[512];
  node_stack[0] = 0;

#ifdef ENABLE_METRICS
  unsigned int node_index_old = 0;
#endif

  // Init isect info as no hit
  intersector.Update(hit_t, static_cast<unsigned int>(-1));

  intersector.PrepareTraversal(ray, options);

  int dir_sign[3];
  dir_sign[0] = ray.dir[0] < static_cast<T>(0.0) ? 1 : 0;
  dir_sign[1] = ray.dir[1] < static_cast<T>(0.0) ? 1 : 0;
  dir_sign[2] = ray.dir[2] < static_cast<T>(0.0) ? 1 : 0;

  real3<T> ray_inv_dir;
  real3<T> ray_dir;
  ray_dir[0] = ray.dir[0];
  ray_dir[1] = ray.dir[1];
  ray_dir[2] = ray.dir[2];

  ray_inv_dir = vsafe_inverse(ray_dir);

  real3<T> ray_org;
  ray_org[0] = ray.org[0];
  ray_org[1] = ray.org[1];
  ray_org[2] = ray.org[2];

  T min_t = std::numeric_limits<T>::max();
  T max_t = -std::numeric_limits<T>::max();

  if (history != NULL)
    {
    history->clear();
    history->reserve(1024);///(10 * sizeof(node_stack) / sizeof(node_stack[0]));  // Likely 10 * (tree depth) would be OK
    }

  while (node_stack_index >= 0) {

#ifdef ENABLE_METRICS
    stats.NC  += 1;
    stats.BLB += sizeof(BVHNode<T>);
#endif

    unsigned int index = node_stack[node_stack_index];
    const BVHNode<T>& node = BVHAccel<T>::nodes_[index];
#ifdef ENABLE_METRICS
    for (int i = 0; i < TREELET_ARR_SIZE; i++)
      {
      if (std::abs(int(index - node_index_old)) * sizeof(BVHNode<T>) >= size_t(treelet_sizes[i]))
        stats.LJC[i]++;
      }
    node_index_old = index;
#endif

    node_stack_index--;
#ifdef ENABLE_METRICS
    stats.SOC++;
    stats.SBL += sizeof(int);
#endif


    // 'hit' receives 0 if ray missed and 1,2,3 if it hit box wall orthogonal to X,Y,Z respectively
    int hit = IntersectRayAABBExt(&min_t, &max_t, ray.min_t, hit_t, node.bmin,
                                  node.bmax, ray_org, ray_inv_dir, dir_sign);

    if (history != NULL && (hit || index != 0)) // if ray missed the root (node 0), do not count it
      {
      int value = (((hit + 1) << 24) & END_MASK) | (index & START_MASK);
      int side = (value & END_MASK) >> 24;
      Assert(side == hit + 1);
      history->emplace_back(value);
      }

    if (hit) {
      // Branch node
      if (node.flag == 0) {
        int order_near = dir_sign[node.axis];
        int order_far = 1 - order_near;

        // Traverse near first.
        node_stack[++node_stack_index] = node.data[order_far];
        node_stack[++node_stack_index] = node.data[order_near];

#ifdef ENABLE_METRICS
        stats.SOC += 2;
        stats.SBL += 2 * sizeof(int);
#endif
        }
      else {
#ifdef ENABLE_METRICS
        unsigned int num_primitives = node.data[0];
        stats.LC++;
        stats.LC2++;
        stats.TC += num_primitives;
        stats.BLB += num_primitives * (3 * sizeof(uint32_t) + 9 * sizeof(float));
#endif

        if (BVHAccel<float>::TestLeafNode(node, ray, intersector)) {  // Leaf node
        hit_t = intersector.GetT();
        }
      }
    }
    }

  assert(node_stack_index < kNANORT_MAX_STACK_DEPTH);

  bool hit = (intersector.GetT() < ray.max_t);
  intersector.PostTraversal(ray, hit, isect);

  return hit;
  }




/// The number of nodes in the subtree
template <typename T>
inline int BVHAccelExt<T>::NumSubTreeNodes(int root_node_idx) const
  {
  const BVHNode<T>& node = BVHAccel<T>::nodes_[root_node_idx];

  int num_nodes = 1; // root
  if (!node.flag)
    num_nodes += NumSubTreeNodes(node.data[0]) + NumSubTreeNodes(node.data[1]);
  return num_nodes;
  }

/// Determines which nodes (of a sub-tree) constitute its root cluster
/// and which (those "below" the cluster in the tree hierarchy) will be 
/// the roots of child clusters. 
/// 
/// Both arrays of nodes are sorted so that the "criterion" w0 (can be arbitrary
/// but aimed mainly for probability of ray hitting the node box) is descending.
/// 
/// @param[in] root_node_idx - index of the root of this subtree
/// @param[in] w0            - weights of nodes ("criterion")
/// @param[in] threshold     - threshold to select nodes for the root cluster
/// @param[in] max_num_nodes - max allowed number of cluster nodes (but it can contain less)
/// @param[out] rcluster_idx - indices of the nodes that will form the root cluster
///                            Sorted according to criterion
/// @param[out] child_cluster_root_idx - indices of the nodes that will become roots of
///                                      the child clusters.
///                                      Sorted according to criterion
/// @return SUCCESS/FAILURE
template <typename T>
inline int BVHAccelExt<T>::RootCluster4Subtree(int root_node_idx,
  const std::vector<double>& w0,
  double threshold,
  int max_num_nodes,
  std::vector<int>& rcluster_idx,
  std::vector<int>& child_cluster_root_idx) const
  {
  std::vector<int> last;
  std::vector<FloatAndInt> w0_next, w0_child_cluster_root;

  rcluster_idx.reserve(max_num_nodes); // Better max. expected number of nodes in cluster
  rcluster_idx.clear();


  int iter;
  // 'last' contains the deepest level nodes ALREADY included in the cluster 

  last.reserve(rcluster_idx.size() / 2); // For binary tree = approx. the count of nodes in the "lower level" of cluster
  last.clear();
  last.emplace_back(root_node_idx);

  double sum = 1e99;
  for (iter = 0; iter < 1000; iter++)
    {
    if (rcluster_idx.capacity() < rcluster_idx.size() + last.size())
      rcluster_idx.reserve(rcluster_idx.size() + last.size());
    for (size_t j = 0; j < last.size(); j++)
      rcluster_idx.emplace_back(last[j]);

    // Cluster includes nodes from 'last'.

    // 'sum' contains the probability that a ray hits ANY of the 'last' nodes.
    // If it is too low, do not include these nodes into root cluster, but
    // not loose them, add to the child clusters
    const bool end = (sum < 0.01 * w0[root_node_idx]);

    // This iteration will try to add to the cluster some (or all) of their child nodes.

    w0_next.reserve(2 * last.size());
    w0_next.clear();
    sum = 0;
    for (size_t i = 0; i < last.size(); i++)
      {
      const BVHNode<T>& node = BVHAccel<T>::nodes_[last[i]];
      if (!node.flag) // it is not a terminal (=leaf) node
        {
        const int& ic1 = node.data[0], & ic2 = node.data[1]; // children
        w0_next.emplace_back(FloatAndInt((float)w0[ic1], ic1));
        w0_next.emplace_back(FloatAndInt((float)w0[ic2], ic2));
        sum += w0[ic1] + w0[ic2];
        }
      }

    // Sort in DESCENDING order in "w0" value
    std::sort(w0_next.begin(), w0_next.end(), FloatAndInt::greater);

    // Leave only those nodes which comprize the "main fraction" of probability
    last.clear();
    const int i_max = max_num_nodes - rcluster_idx.size() - 1;
    double sum1 = 0;
    for (size_t i = 0; i < w0_next.size(); i++)
      {
      const FloatAndInt& w0n = w0_next[i];

      if (!end && int(i) < i_max && sum1 <= threshold * sum) // Remained nodes are not important...
        last.emplace_back(w0n.i); // w0n.i is the current index of node (in their current array)
      else // the node is outside cluster => it is the root of an out-of-cluster subtree
        w0_child_cluster_root.emplace_back(w0n); // weight and index of this child root
      sum1 += w0n.x;
      }
    if (end)
      break;
    } // end for (iter = 0; iter < 1000; iter++)

    // Sort roots of child clusters according to their probability (=w0)
  std::sort(w0_child_cluster_root.begin(), w0_child_cluster_root.end(), FloatAndInt::greater);

  child_cluster_root_idx.resize(w0_child_cluster_root.size());
  for (size_t i = 0; i < w0_child_cluster_root.size(); i++)
    child_cluster_root_idx[i] = w0_child_cluster_root[i].i;

  // Array 'rcluster_idx' is already duly ordered, so:
  return SUCCESS;
  }


/// Suggests ordering of the nodes of a subtree "by clusters". First comes the
/// "root cluster" then "child clusters". Each of the latter, in turn, consists
/// of its root cluster, past which come child clusters, of which... etc.
/// 
/// If a cluster is NOT further subdivided, its nodes are ordered so that their
/// criterion "w0" is descending. 
/// If it is further clusterized, then first comes its root cluster, then child 
/// clusters in such an order that first comes the one whose root has the largest 
/// "w0", then the one whose root has the next "w0" and so on. 
///
/// The order of nodes inside each sub-cluster (root or child) obeys the same
/// obey the same rule, i.e. if this sub-cluster is further subdivided, then... else
/// ordered so that their "w0" is descending.
/// 
/// @param[in] root_node_idx - index of the root of this subtree
/// @param[in] w0            - weights of nodes ("criterion")
/// @param[in] threshold     - threshold to select nodes for the root cluster
/// @param[in] max_num_nodes - max allowed number of cluster nodes (but it can contain less)
/// @param[out] indices      - indices of ALL nodes of the sub-tree in the new order
/// 
/// @return SUCCESS/FAILURE
template <typename T>
inline int BVHAccelExt<T>::ClusterizeSubtree(int root_node_idx,
  const std::vector<double>& w0,
  double threshold,
  std::vector<int>& indices) const
  {
  std::vector<int> child_cluster_root_idx, next;
  int max_num_nodes;

  // Create root cluster
  max_num_nodes = std::min(128, (int)sqrt(double(NumSubTreeNodes(root_node_idx))));///????;
///  max_num_nodes = 10;
  if (RootCluster4Subtree(root_node_idx, w0, threshold, max_num_nodes, indices, child_cluster_root_idx) != SUCCESS)
    {
    Assert(false);
    return FAILURE;
    }

  
  // Create child clusters. Notice the order of child clusters is dictated by
  // 'child_cluster_root_idx' which is sorted according to the "cost criterion"
  for (size_t i = 0; i < child_cluster_root_idx.size(); i++)
  {
    if (ClusterizeSubtree(child_cluster_root_idx[i], w0, threshold, next) != SUCCESS)
      {
      Assert(false);
      return FAILURE;
      }
    if (indices.capacity() < indices.size() + next.size())
      indices.reserve(indices.size() + next.size());
    for (size_t j = 0; j < next.size(); j++)
      indices.emplace_back(next[j]);
  }

  /// Test
  const int num_nodes = NumSubTreeNodes(root_node_idx);
  if (num_nodes != int(indices.size()))
  {
    Assert(false);
    return FAILURE;
  }
  return SUCCESS;
}

/// Changes the order of the tree nodes in memory WHLE RETAINING THE TREE TOPOLOGY
/// and updates the links to the children kept in the nodes accordingly.
/// 
/// New order of nodes is that they are stored "cluster by cluster", each cluster,
/// in turn, is usually further subdivided by sub-clusters etc. But if it is too small
/// it is not further subdivided.
/// 
/// If a cluster is NOT further subdivided, its nodes are ordered so that their
/// criterion "w0" is descending. 
/// If it is further clusterized, then first comes its root cluster, then child 
/// clusters in such an order that first comes the one whose root has the largest 
/// "w0", then the one whose root has the next "w0" and so on. 
///
/// The order of nodes inside each sub-cluster (root or child) obeys the same
/// obey the same rule, i.e. if this sub-cluster is further subdivided, then... else
/// ordered so that their "w0" is descending.  
/// 
/// @param[in] w0            - weights of nodes ("criterion")
/// @param[in] threshold     - threshold to select nodes for the root cluster
/// @return SUCCESS/FAILURE
template <typename T>
inline int BVHAccelExt<T>::ReArrangeByClusters(const std::vector<double>& w0, double threshold)
  {
  std::vector<int> old2new_indices, new2old_indices;

  // New order of the distribution of tree nodes in memory. 
  // The current tree/nodes data itself is not YET changed
  if (ClusterizeSubtree(0, w0, threshold, new2old_indices) != SUCCESS)
    {
    Assert(false);
    return FAILURE;
    }
  Assert(new2old_indices.size() == BVHAccel<T>::nodes_.size());
  old2new_indices.resize(new2old_indices.size());

  std::vector<BVHNode<T> > new_nodes;
  new_nodes.resize(BVHAccel<T>::nodes_.size());

#if 1
  ///  For test
  FILE* fd = fopen("new2old.indices", "w");
  for (size_t i = 0; i < BVHAccel<T>::nodes_.size(); i++)
    fprintf(fd, "%15i %15i %15.5f\n", int(i), new2old_indices[i], w0[i]);
  fclose(fd);
#endif

  // Array of nodes in new order:
  for (size_t i = 0; i < BVHAccel<T>::nodes_.size(); i++)
    {
    const int idx = new2old_indices[i];
    old2new_indices[idx] = i;
    new_nodes[i] = BVHAccel<T>::nodes_[idx];
    }

  BVHAccel<T>::nodes_ = new_nodes;

  // Update links to children so that they fit new array
  for (size_t i = 0; i < BVHAccel<T>::nodes_.size(); i++)
    {
    BVHNode<T>& node = BVHAccel<T>::nodes_[i];
    if (!node.flag) // it is not a terminal (=leaf) node
      {
      node.data[0] = old2new_indices[node.data[0]];
      node.data[1] = old2new_indices[node.data[1]];
      }
    }

  return SUCCESS;
  }

/// Traverse subtree and for small-size nodes calculate frequency of hits by scaling that of parent
///
/// For all nodes of the subtree takes checks the number of ray hits and if too low,
/// calculate it by scaling the number of hits for the parent node by the ratio
/// of areas.
/// 
/// @param[in]     root_node_idx - index of the root of this subtree
/// @param[in,out] num_hits      - number of hits, separately for the 3 walls normal to X,Y,Z
template <typename T>
inline int BVHAccelExt<T>::ImproveSmallNodesProb(int root_node_idx, std::vector< real3<double> >& num_hits) const
  {
  static const int low_hit_thereshold = 500; // For more hits retain raw data

  const BVHNode<T>& node = BVHAccel<T>::nodes_[root_node_idx];
  if (node.flag) // terminal node, no children
    return SUCCESS;

  // Child nodes:
  const int ic1 = node.data[0], ic2 = node.data[1]; 
  const BVHNode<T> &child1 = BVHAccel<T>::nodes_[ic1], 
                   &child2 = BVHAccel<T>::nodes_[ic2];

  // Number of ray hits by the node's box walls (separately for walls normal to X, Y and Z)
  const real3<double>& n_root = num_hits[root_node_idx];
  real3<double> &nc1 = num_hits[ic1], &nc2 = num_hits[ic2];

  const real3<double> &s = node.Areas();

  if (nc1[0] + nc1[1] + nc1[2] < low_hit_thereshold) // else retain raw data
    {
    const real3<double>& s1 = child1.Areas();
    // correct child:
    nc1[0] = n_root[0] * s1[0] / s[0];
    nc1[1] = n_root[1] * s1[1] / s[1];
    nc1[2] = n_root[2] * s1[2] / s[2];
    }
  // and its subtree:
  if (ImproveSmallNodesProb(ic1, num_hits) != SUCCESS)
    {
    Assert(false);
    return FAILURE;
    }

  if (nc2[0] + nc2[1] + nc2[2] < low_hit_thereshold) // else retain raw data
    {
    const real3<double>& s2 = child2.Areas();
    // correct child:
    nc2[0] = n_root[0] * s2[0] / s[0];
    nc2[1] = n_root[1] * s2[1] / s[1];
    nc2[2] = n_root[2] * s2[2] / s[2];
    }
  // and its subtree:
  if (ImproveSmallNodesProb(ic2, num_hits) != SUCCESS)
    {
    Assert(false);
    return FAILURE;
    }

  return SUCCESS;
  }


#ifdef __clang__
#pragma clang diagnostic pop
#endif
}  // namespace nanortext

#endif  // NANORTEXT_H_
