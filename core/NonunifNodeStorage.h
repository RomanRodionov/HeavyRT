#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>

#ifdef HALFFLOAT
#include "../external/halffloat/half2.2/half.hpp"
#endif

#include "LiteMath.h"
#include "builders/FatBVH.h"

using LiteMath::BBox3f;
using LiteMath::float3;
using LiteMath::int4;
using LiteMath::uint;

constexpr uint MAX_OFFSET_MASK     = 0x07FFFFFF;  // max 27 bits for offset
constexpr uint MAX_TRG_START_MASK  = 0x007FFFFF;  // max 23 bits for the first triangle index
constexpr uint MAX_TRG_COUNT_MASK  = 0x0000000F;  // max 4 bits for the number of triangles in a leaf

#ifdef HALFFLOAT
class BoxScalerHalf
  {
  public:
    /// Converts the "standard/exact" input into internal form
    void Set(const LiteMath::float3 &origin, const LiteMath::float3& size);

    /// Converts the internal form into "standard/exact" values
    void Get(LiteMath::float3& origin, LiteMath::float3& ssize) const;

    float GetOutSizeScale() const { return 63.0f; }

  protected:
    /// Box origin; the "w" component is not used and is added for 16bytes alignment
    half_float::half m_origin_x, m_origin_y, m_origin_z, m_origin_w;
    /// Box size/63.0 where 63 is the range for 6 bit numbers; the "w" component is not used and is added for 16bytes alignment
    half_float::half m_ssize_x, m_ssize_y, m_ssize_z, m_ssize_w;
  };
#endif

class NonunifNodeStorage
  {
  public:
    /// How many data elements (of size of the own node data) the treelet's bounding box takes
    inline static int TreeletBBoxTakesElements();
  public:
    /// Intialize from a "full" (not compressed) treeletized/clusterized tree.
    int Create(const std::vector<BVHNodeFat>& a_nodes, const std::vector<uint32_t> &a_bvhOffsets, const std::vector<int>& a_treelet_root, const std::vector<int>& grStart);

    /// Appends the data, e.g. for a new geometry objected added to the scene
    int Append(const NonunifNodeStorage &nns);

    /// Given the pointer to the compressed node data, calculates the bounding boxes of its left and right children
    /// The input data can relate to the "internal" or "root" node.
    /// It takes the pointer to the treelet root, shifts it left by 2*sizeof(float4) and
    /// this is the pointer to the "pivot" bounding box stored as two float4.
    /// Then it extracts from the current node 2x3 6bit numbers that represent
    /// the left bounding box, RELATIVE to that of the root. Same for the right
    /// bounding box. Scales then by the "pivot" box, and voila
    void GetChildBoxes(const LiteMath::uint4* data, BBox3f& left, BBox3f& right) const;

    /// Another variant, where the "pivot" bounding box (to scale the 6bit relative data)
    /// is not found internally, but is provided explcitly as an argument.
    static void GetChildBoxes(const LiteMath::uint4* data, const LiteMath::float4 bb_data[2], BBox3f& left, BBox3f& right);

    /// Again the "pivot" bounding box (to scale the 6bit relative data)
    /// is not found internally, but is provided explcitly as an argument, but as float3.
    static void GetChildBoxes(const LiteMath::uint4* data, const LiteMath::float3 &bb_origin, const LiteMath::float3& bb_ssize, BBox3f& left, BBox3f& right);

    /// Gets offs_left/offs_right in 32 bit form.
    static void GetOffsets(uint& offs_left, uint& offs_right, const LiteMath::uint4& data);

    /// Extracts the starting position of the treelet root
    inline const LiteMath::uint4* GetTreeletRoot(const LiteMath::uint4* data) const;

    /// Extracts the data on the bounding box of the treelet's root
    inline static void GetRootBBox(LiteMath::float3 &origin, LiteMath::float3& ssize, const LiteMath::uint4* tr_root);

    /// Sets the data on the bounding box of the treelet's root
    inline static void SetRootBBox(const LiteMath::BBox3f& bbox, LiteMath::uint4* tr_root);

    /// Offset of the k-th BVH tree in the global data
    inline uint32_t Offset(int k) const
      {
      return m_bvhOffsets16B[k];
      }

    /// Gets the data of the node which starts from current position 'pos'
    inline const LiteMath::uint4* operator[](uint pos) const
      {
      return (data.data() + pos);
      }

    size_t DataSizeInBytes() const { return data.size()*sizeof(uint4); }

  protected:
    /// @name Basic functions to nodes concatenated in a uint4 array
    //@{
    /// Extracts from the current node 2x3 6bit numbers that represent the left
    /// bounding box (RELATIVE to that of the root), and 2x3 numbers for the right box.
    /// These numbers are then put in arrays of bytes: first 3 values for (x,y,z) of 
    /// the min point, then 3 values for (x,y,z) of the max point of each box.
    static void GetChildBoxes(const LiteMath::uint4* data, LiteMath::uchar left[6], LiteMath::uchar right[6]);
    //@}

    /// @name For creation from "full" nodes
    //@{
    //////////////////////////////////////////////////////////////////////////////
    /// Sets (relative) child bounding boxes in 6bits format
    inline static void SetChildBoxes(const LiteMath::uchar left[6], const LiteMath::uchar right[6], LiteMath::uint4* data);
    
    /// Sets (relative) child bounding boxes in 6bits format
    inline static void SetChildBoxes(const LiteMath::BBox3f &cb1, const LiteMath::BBox3f &cb2, const LiteMath::float3 &tr_bb_origin, const LiteMath::float3& tr_bb_rev_ssize, LiteMath::uint4* data);

    //////////////////////////////////////////////////////////////////////////////
    /// Sets offs_left/offs_right in 25 bit form. The input values must include flags "leaf"
    static void SetOffsets(uint offs_left, uint offs_right, LiteMath::uint4* data);
    //@}

  public:
    /// Data (bounding boxes for the roots of treelets and the nodes themself)
    /// Can merge several trees (for several geometric objects)
    std::vector<LiteMath::uint4> data;
    /// Positions of the individual BVH trees in 'data'.
    /// The data for the k-th tree starts at (*this)[a_bvhOffsets16B[k]]. 
    /// First come 2xfloat4 (exactly at this position) that represent the bounding box 
    /// of the tree root, then the root node itself, etc.
    std::vector<uint32_t> m_bvhOffsets16B;

    /// log2(1 + number of nodes in group)
    int m_log2_group_size;
  };

/// How many data elements (of size of the own node data) the treelet's bounding box takes
int NonunifNodeStorage::TreeletBBoxTakesElements()
  {
#ifdef HALFFLOAT
  // BoxScalerHalf = 8 halffloat values
  return 1;
#else
  // Two float4 values
  return 2;
#endif
  }

/// Extracts the starting position of the treelet root
/// It is calculated by shifting "downwards" the pointer to current data 'data'
/// and returns pointer to the base data of the root of current treelet.
/// The root's extra data (now two float4 for bounding box) PRECEEDS it (has smaller pointer).
/// @param[in] data  - node data
const LiteMath::uint4* NonunifNodeStorage::GetTreeletRoot(const LiteMath::uint4* a_data) const
  {
//  return 1 + this->data.data() + 8 * ((data - this->data.data()) / 8);
  return 1 + this->data.data() + (((a_data - this->data.data()) >> m_log2_group_size) << m_log2_group_size);
  }

/// Extracts the data on the bounding box of the treelet's root
/// @param[out] origin - min point of the bounding box of the treelet's root
/// @param[out] ssize  - size/63 of the bounding box of the treelet's root
/// @param[in] tr_root  - root node of treelet
void NonunifNodeStorage::GetRootBBox(LiteMath::float3& origin, LiteMath::float3& ssize, const LiteMath::uint4* tr_root)
  {
#ifdef HALFFLOAT
  const BoxScalerHalf* bbscaler = reinterpret_cast<const BoxScalerHalf*>(tr_root) - 1; // '1'=sizeof(BoxScalerHalf)/sizeof(uint4)
  bbscaler->Get(origin, ssize);
#else
  static const float scale63 = (1.0 / 63.0); // '63' is teh range of 6bit numbers

  // Pointer to the treelet's root bounding box origin
  const LiteMath::float4* tr_bbox_orig = reinterpret_cast<const LiteMath::float4*>(tr_root) - 2; // '2'=(2*sizeof(float4))/sizeof(uint4)

  // Pointer to the bounding box scaled size
  const LiteMath::float4* tr_bbox_ssize = tr_bbox_orig + 1;

  origin.x = tr_bbox_orig->x;
  origin.y = tr_bbox_orig->y;
  origin.z = tr_bbox_orig->z;

  ssize.x = tr_bbox_ssize->x;
  ssize.y = tr_bbox_ssize->y;
  ssize.z = tr_bbox_ssize->z;
#endif
  }

/// Sets the data on the bounding box of the treelet's root
/// @param[in] bbox    - bounding box (comprisizing whole node, i.e. left and right child boxes!)
/// @param[in] tr_root - root node of treelet
void NonunifNodeStorage::SetRootBBox(const LiteMath::BBox3f &bbox, LiteMath::uint4* tr_root)
  {
#ifdef HALFFLOAT
  BoxScalerHalf* bbscaler = reinterpret_cast<BoxScalerHalf*>(tr_root) - 1; // '1'=sizeof(BoxScalerHalf)/sizeof(uint4)
  bbscaler->Set(bbox.boxMin, bbox.boxMax - bbox.boxMin);
#else
  static const float scale63 = (1.0 / 63.0); // '63' is teh range of 6bit numbers

  // Pointer to the treelet's root bounding box origin
  LiteMath::float4* tr_bbox_orig = reinterpret_cast<LiteMath::float4*>(tr_root) - 2; // '2'=(2*sizeof(float4))/sizeof(uint4)

  // Pointer to the bounding box scaled size
  LiteMath::float4* tr_bbox_ssize = tr_bbox_orig + 1;

  tr_bbox_orig->x = bbox.boxMin.x;
  tr_bbox_orig->y = bbox.boxMin.y;
  tr_bbox_orig->z = bbox.boxMin.z;

  tr_bbox_ssize->x = (bbox.boxMax.x - bbox.boxMin.x) * scale63;
  tr_bbox_ssize->y = (bbox.boxMax.y - bbox.boxMin.y) * scale63;
  tr_bbox_ssize->z = (bbox.boxMax.z - bbox.boxMin.z) * scale63;
#endif
  }


//////////////////////////////////////////////////////////////////////////////
/// Sets (relative) child bounding boxes in 6bits format
void NonunifNodeStorage::SetChildBoxes(const LiteMath::uchar left[6], const LiteMath::uchar right[6], LiteMath::uint4* data)
  {
  data->M[0] &= 0x00000000;  // zeroing data[0]: bits:[0:31]
  data->M[1] &= 0x00000000;  // zeroing data[1]: bits:[0:31]
  data->M[2] &= 0xFFFFFF00;  // zeroing data[2]: bits:[0:7]

  for (int i = 0; i < 6; i++)
    {
    data->M[0] |= ((uint(left[i])) << (i * 6));
    }
  data->M[1] |= ((uint(left[5])) >> 2);  // put remaining four bits of lmax_z into data[1]: bits:[0:3]

  for (int i = 0; i < 5; i++)
    {
    data->M[1] |= ((uint(right[i])) << (i * 6 + 4));
    }
  data->M[2] |= ((uint(right[4])) >> 4);  // put remaining two bits of rmax_y into data[2]: bits:[0:1]
  data->M[2] |= ((uint(right[5])) << 2);  // also put rmax_z into data[2]: bits:[2:7]
  }

//////////////////////////////////////////////////////////////////////////////
/// Sets (relative) child bounding boxes in 6bits format
/// 
/// @param[in] cb1 - bounding box of the left child
/// @param[in] cb2 - bounding box of the right child
/// @param[out] tr_bb_origin     - min point of the bounding box of the treelet's root
/// @param[out] tr_bb_rev_ssize  - 63/size of the bounding box of the treelet's root
/// @param[in,out] data
void NonunifNodeStorage::SetChildBoxes(const LiteMath::BBox3f& cb1, const LiteMath::BBox3f& cb2, const LiteMath::float3& tr_bb_origin, const LiteMath::float3& tr_bb_rev_ssize, LiteMath::uint4* data)
  {
  // Calculate RELATIVE (to B.B. of the treelet root) bounding boxes and convert to
  // 6bits format
  LiteMath::uchar left[6], right[6];
  left[0] = (std::min(63.0f, std::max(0.0f, (cb1.boxMin[0] - tr_bb_origin[0]) * tr_bb_rev_ssize[0])));
  left[1] = (std::min(63.0f, std::max(0.0f, (cb1.boxMin[1] - tr_bb_origin[1]) * tr_bb_rev_ssize[1])));
  left[2] = (std::min(63.0f, std::max(0.0f, (cb1.boxMin[2] - tr_bb_origin[2]) * tr_bb_rev_ssize[2])));
  left[3] = (std::min(63.0f, std::max(0.0f, std::ceil((cb1.boxMax[0] - tr_bb_origin[0]) * tr_bb_rev_ssize[0]))));
  left[4] = (std::min(63.0f, std::max(0.0f, std::ceil((cb1.boxMax[1] - tr_bb_origin[1]) * tr_bb_rev_ssize[1]))));
  left[5] = (std::min(63.0f, std::max(0.0f, std::ceil((cb1.boxMax[2] - tr_bb_origin[2]) * tr_bb_rev_ssize[2]))));

  right[0] = (std::min(63.0f, std::max(0.0f, (cb2.boxMin[0] - tr_bb_origin[0]) * tr_bb_rev_ssize[0])));
  right[1] = (std::min(63.0f, std::max(0.0f, (cb2.boxMin[1] - tr_bb_origin[1]) * tr_bb_rev_ssize[1])));
  right[2] = (std::min(63.0f, std::max(0.0f, (cb2.boxMin[2] - tr_bb_origin[2]) * tr_bb_rev_ssize[2])));
  right[3] = (std::min(63.0f, std::max(0.0f, std::ceil((cb2.boxMax[0] - tr_bb_origin[0]) * tr_bb_rev_ssize[0]))));
  right[4] = (std::min(63.0f, std::max(0.0f, std::ceil((cb2.boxMax[1] - tr_bb_origin[1]) * tr_bb_rev_ssize[1]))));
  right[5] = (std::min(63.0f, std::max(0.0f, std::ceil((cb2.boxMax[2] - tr_bb_origin[2]) * tr_bb_rev_ssize[2]))));

  SetChildBoxes(left, right, data);

  //#define OSE_TESTS1
  //#ifdef OSE_TESTS1 // tests
  //  assert(left[0] >= 0 && left[0] <= 63);
  //  assert(left[1] >= 0 && left[1] <= 63);
  //  assert(left[2] >= 0 && left[2] <= 63);
  //  assert(left[3] >= 0 && left[3] <= 63);
  //  assert(left[4] >= 0 && left[4] <= 63);
  //  assert(left[5] >= 0 && left[5] <= 63);
  //
  //  assert(right[0] >= 0 && right[0] <= 63);
  //  assert(right[1] >= 0 && right[1] <= 63);
  //  assert(right[2] >= 0 && right[2] <= 63);
  //  assert(right[3] >= 0 && right[3] <= 63);
  //  assert(right[4] >= 0 && right[4] <= 63);
  //  assert(right[5] >= 0 && right[5] <= 63);
  //
  //  LiteMath::BBox3f cb1_, cb2_;
  //  GetChildBoxes(data, cb1_, cb2_);
  //  static const double eps = 1e-5;
  //  assert(cb1_.boxMin.x <= cb1.boxMin.x + eps && cb1_.boxMax.x >= cb1.boxMax.x - eps);
  //  assert(cb1_.boxMin.y <= cb1.boxMin.y + eps && cb1_.boxMax.y >= cb1.boxMax.y - eps);
  //  assert(cb1_.boxMin.z <= cb1.boxMin.z + eps && cb1_.boxMax.z >= cb1.boxMax.z - eps);
  //  assert(cb2_.boxMin.x <= cb2.boxMin.x + eps && cb2_.boxMax.x >= cb2.boxMax.x - eps);
  //  assert(cb2_.boxMin.y <= cb2.boxMin.y + eps && cb2_.boxMax.y >= cb2.boxMax.y - eps);
  //  assert(cb2_.boxMin.z <= cb2.boxMin.z + eps && cb2_.boxMax.z >= cb2.boxMax.z - eps);
  //#endif
  }
