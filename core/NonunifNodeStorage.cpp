#include <cassert>
#include "NonunifNodeStorage.h"

#define OSE_TESTS // For self-tests. Decelerates calculations


/// Calculates common bounding box for the given group of
/// nodes. It is the union of BBs of all that nodes.
/// @param[in] a_nodes - array of nodes, may merge several BVH tree
/// @param[in] from    - the index of the first node in the group in 'a_nodes'
/// @param[in] to      - the index of the last  node in the group in 'a_nodes'
/// @return Resulting box
static LiteMath::BBox3f CalcCommonBB(const std::vector<BVHNodeFat>& a_nodes,
                                     uint32_t from, uint32_t to)
  {
  LiteMath::BBox3f bbox;
  bbox.boxMin.x =  1e30;
  bbox.boxMin.y =  1e30;
  bbox.boxMin.z =  1e30;
  bbox.boxMax.x = -1e30;
  bbox.boxMax.y = -1e30;
  bbox.boxMax.z = -1e30;

  for (uint32_t i = from; i <= to; i++)
    {
    const BVHNodeFat& node = a_nodes[i];

    const BBox3f cb1 = GetChildBoxLeft(node);
    const BBox3f cb2 = GetChildBoxRight(node);

    bbox.boxMin.x = std::min(bbox.boxMin.x, std::min(cb1.boxMin.x, cb2.boxMin.x));
    bbox.boxMin.y = std::min(bbox.boxMin.y, std::min(cb1.boxMin.y, cb2.boxMin.y));
    bbox.boxMin.z = std::min(bbox.boxMin.z, std::min(cb1.boxMin.z, cb2.boxMin.z));
    bbox.boxMax.x = std::max(bbox.boxMax.x, std::max(cb1.boxMax.x, cb2.boxMax.x));
    bbox.boxMax.y = std::max(bbox.boxMax.y, std::max(cb1.boxMax.y, cb2.boxMax.y));
    bbox.boxMax.z = std::max(bbox.boxMax.z, std::max(cb1.boxMax.z, cb2.boxMax.z));
    }
  return bbox;
  }

#ifdef HALFFLOAT

static inline half_float::half roundToNeg(float a_val) { return half_float::nextafter(half_float::half(a_val), -std::numeric_limits<half_float::half>::max()); }
static inline half_float::half roundToPos(float a_val) { return half_float::nextafter(half_float::half(a_val), +std::numeric_limits<half_float::half>::max()); }

/// Converts the "standard/exact" input into internal form
/// 
/// Ensures that these rounded to half-float values are rounded so that
/// the related bounding box @b includes the original one (i.e. is not smaller)
/// But it can be larger.
/// 
/// @param[in] origin - the "minimal point" of bounding box
/// @param[in] size   - size of the bounding box
void BoxScalerHalf::Set(const LiteMath::float3& origin, const LiteMath::float3& size)
{
  static const float scale63 = (1.0 / 63.0); // '63' is teh range of 6bit numbers

  // We must ensure that the low-precision box INCLUDES the exact one. For this purpose,
  // its origin must be rounded towards the nearest LOWER value
  //m_origin_x = half_float::half_cast<half_float::half, std::round_toward_neg_infinity, float>(origin.x);
  //m_origin_y = half_float::half_cast<half_float::half, std::round_toward_neg_infinity, float>(origin.y);
  //m_origin_z = half_float::half_cast<half_float::half, std::round_toward_neg_infinity, float>(origin.z);
  m_origin_x = roundToNeg(origin.x);
  m_origin_y = roundToNeg(origin.y);
  m_origin_z = roundToNeg(origin.z);

  // With the box size, situation is a bit more complex. We must ensure that
  // origin2+size2 >= origin+size, where origin2 and size2 are the
  // 2byte representations. 
  // origin2 has been defined. It is SMALLER than origin. Thus even if we round 
  // ssize UP, it could happen that origin2+size2 < origin+size
  // Thus we must first re-calculate the exact value of ssize and only then
  // convert it.
  const LiteMath::float3 bboxMax(origin + size);
  const LiteMath::float3 origin2(m_origin_x, m_origin_y, m_origin_z);
  const LiteMath::float3 ssize2 = (bboxMax - origin2) * scale63;

  //m_ssize_x = half_float::half_cast<half_float::half, std::round_toward_infinity, float>(ssize2.x);
  //m_ssize_y = half_float::half_cast<half_float::half, std::round_toward_infinity, float>(ssize2.y);
  //m_ssize_z = half_float::half_cast<half_float::half, std::round_toward_infinity, float>(ssize2.z);
  m_ssize_x = roundToPos(ssize2.x);
  m_ssize_y = roundToPos(ssize2.y);
  m_ssize_z = roundToPos(ssize2.z);
}

/// Converts the internal form into "standard/exact" values
/// @param[in] origin - the "minimal point" of bounding box
/// @param[in] ssize  - @b scaled size of the bounding box (scale/63). The scaling 
///                     is needed because these values are intended for multiplication
///                     by 6bit (range = 63) encodings of values in [0,1].
void BoxScalerHalf::Get(LiteMath::float3& origin, LiteMath::float3& ssize) const
  {
  origin.x = m_origin_x;
  origin.y = m_origin_y;
  origin.z = m_origin_z;

  ssize.x = m_ssize_x;
  ssize.y = m_ssize_y;
  ssize.z = m_ssize_z;
  }
#endif

/// Extracts from the current node 2x3 6bit numbers that represent the left
/// bounding box (RELATIVE to that of the root), and 2x3 numbers for the right box.
/// These numbers are then put in arrays of bytes: first 3 values for (x,y,z) of 
/// the min point, then 3 values for (x,y,z) of the max point of each box.
void NonunifNodeStorage::GetChildBoxes(const LiteMath::uint4* data,
                                       LiteMath::uchar left[6], LiteMath::uchar right[6])
  {
  /*
  constexpr uint LMINX_MASK   = 0x0000003F; data[0]: bits:[0:5]
  constexpr uint LMINY_MASK   = 0x00000FC0; data[0]: bits:[6:11]
  constexpr uint LMINZ_MASK   = 0x0003F000; data[0]: bits:[12:17]
  constexpr uint LMAXX_MASK   = 0x00FC0000; data[0]: bits:[18:23]
  constexpr uint LMAXY_MASK   = 0x3F000000; data[0]: bits:[24:29]
  constexpr uint LMAXZ_MASK_L = 0xC0000000; data[0]: bits:[30:31] 2 bit
  constexpr uint LMAXZ_MASK_H = 0x0000000F; data[1]: bits:[0:3]   4 bit

  constexpr uint RMINX_MASK   = 0x000003F0; data[1]: bits:[4:9]
  constexpr uint RMINY_MASK   = 0x0000FC00; data[1]: bits:[10:15]
  constexpr uint RMINZ_MASK   = 0x003F0000; data[1]: bits:[16:21]
  constexpr uint RMAXX_MASK   = 0x0FC00000; data[1]: bits:[22:27]
  constexpr uint RMAXY_MASK_L = 0xF0000000; data[1]: bits:[28:31] 4 bit
  constexpr uint RMAXY_MASK_H = 0x00000003; data[2]: bits:[0:1]   2 bit
  constexpr uint RMAXZ_MASK   = 0x000000FC; data[2]: bits:[2:7]
  */

  uint l_mask = 0x0000003F;
  for (int i = 0; i < 6; i++)
    {
    left[i] = (data->M[0] & l_mask) >> (i * 6);
    l_mask = l_mask << 6;
    }
  left[5] = ((uint)left[5]) | (data->M[1] & 0x0000000F) << 2;  // combine from two adjacent ints

  uint r_mask = 0x000003F0;
  for (int i = 0; i < 5; i++)
    {
    right[i] = (data->M[1] & r_mask) >> (i * 6 + 4);
    r_mask = r_mask << 6;
    }
  right[4] = ((uint)right[4]) | (data->M[2] & 0x00000003) << 4;  // combine from two adjacent ints
  right[5] = (data->M[2] & 0x000000FC) >> 2;
  }


//////////////////////////////////////////////////////////////////////////////
/// Sets offs_left/offs_right in 28 bit form. The input values must include flags "leaf"
/// @param[in] offs_left - either the "leaf data" (the same as for BVHNodeFat) or
///                        the offset of the left child node from the tree start
///                        (given by @c Offset(int tree_index)).
///                        That is, this child node is (*this)[offs_left + Offset(tree_index)]
/// @param[in] offs_right - either the "leaf data" (the same as for BVHNodeFat) or
///                         the offset of the right child node from the tree start
///                        (given by @c Offset(int tree_index))
///                        That is, this child node is (*this)[offs_right + Offset(tree_index)]
/// @param[in,out] data  - node data
void NonunifNodeStorage::SetOffsets(uint offs_left, uint offs_right, LiteMath::uint4* data)
  {
  // left part
  uint comp_offset_left    = EXTRACT_OFFSET(offs_left) & MAX_OFFSET_MASK;
  uint comp_trg_start_left = EXTRACT_START(offs_left) & MAX_TRG_START_MASK;
  uint comp_trg_count_left = EXTRACT_COUNT(offs_left) & MAX_TRG_COUNT_MASK;

  data->M[2] &= 0x000000FF;  // zeroing data[2]: bits:[8:31]
  data->M[3] &= 0xFFFFFFF8;  // zeroing data[3]: bits:[0:2]
  if (offs_left & LEAF_BIT)
    {
    assert(!(EXTRACT_START(offs_left) & (~MAX_TRG_START_MASK)));
    assert(!(EXTRACT_COUNT(offs_left) & (~MAX_TRG_COUNT_MASK)));
    data->M[2] |= (comp_trg_start_left << 8);   // data[2]: bits:[8:30]
    data->M[2] |= (comp_trg_count_left << 31);  // data[2]: bits:[31:31]
    data->M[3] |= (comp_trg_count_left >> 1);   // data[3]: bits:[0:2]
    data->M[3] |= 0x00000008;  // set leaf bit (data[3]: bit:3)
    }
  else
    {
    assert(!(EXTRACT_OFFSET(offs_left) & (~MAX_OFFSET_MASK)));
    data->M[2] |= (comp_offset_left << 8);   // data[2]: bits:[8:31]
    data->M[3] |= (comp_offset_left >> 24);  // data[3]: bits:[0:2]
    data->M[3] &= 0xFFFFFFF7;  // clear leaf bit (data[3]: bit:3)
    }

  // right part
  uint comp_offset_right = EXTRACT_OFFSET(offs_right) & MAX_OFFSET_MASK;
  uint comp_trg_start_right = EXTRACT_START(offs_right) & MAX_TRG_START_MASK;
  uint comp_trg_count_right = EXTRACT_COUNT(offs_right) & MAX_TRG_COUNT_MASK;

  data->M[3] &= 0x8000000F;  // zeroing data[3]: bits:[4:30]
  if (offs_right & LEAF_BIT)
    {
    assert(!(EXTRACT_START(offs_right) & (~MAX_TRG_START_MASK)));
    assert(!(EXTRACT_COUNT(offs_right) & (~MAX_TRG_COUNT_MASK)));
    data->M[3] |= (comp_trg_start_right << 4);   // data[3]: bits:[4:26]
    data->M[3] |= (comp_trg_count_right << 27);  // data[3]: bits:[27:30]
    data->M[3] |= 0x80000000;  // set leaf bit (data[3]: bit:31)
    }
  else
    {
    assert(!(EXTRACT_OFFSET(offs_right) & (~MAX_OFFSET_MASK)));
    data->M[3] |= (comp_offset_right << 4);  // data[3]: bits:[4:30]
    data->M[3] &= 0x7FFFFFFF;  // clear leaf bit (data[3]: bit:31)
    }
  }

//////////////////////////////////////////////////////////////////////////////
/// Gets offs_left/offs_right in 32 bit form.
/// @param[in] offs_left - either the "leaf data" (the same as for BVHNodeFat) or
///                        the offset of the left child node from the tree start
///                        (given by @c Offset(int tree_index))
///                        That is, this child node is (*this)[offs_left + Offset(tree_index)]
/// @param[in] offs_right - either the "leaf data" (the same as for BVHNodeFat) or
///                         the offset of the right child node from the tree start
///                        (given by @c Offset(int tree_index))
///                        That is, this child node is (*this)[offs_right + Offset(tree_index)]
/// @param[in,out] data  - node data
void NonunifNodeStorage::GetOffsets(uint& offs_left, uint& offs_right, const LiteMath::uint4& data)
  {
  offs_left = offs_right = 0;
  // left part
  if (data.M[3] & 0x00000008)  // check leaf bit (data[3]: bit:3)
    {
    uint start = 0, count = 0;

    start |= data.M[2] >> 8;
    start &= MAX_TRG_START_MASK;

    count |= data.M[2] >> 31;
    count |= data.M[3] << 1;
    count &= MAX_TRG_COUNT_MASK;

    offs_left = PackOffsetAndSize(start, count);
    }
  else
    {
    offs_left |= data.M[2] >> 8;
    offs_left |= data.M[3] << 24;
    offs_left &= MAX_OFFSET_MASK;
    }

  // right part
  if (data.M[3] & 0x80000000)  // check leaf bit (data[3]: bit:31)
    {
    uint start = 0, count = 0;

    start |= data.M[3] >> 4;
    start &= MAX_TRG_START_MASK;

    count |= data.M[3] >> 27;
    count &= MAX_TRG_COUNT_MASK;

    offs_right = PackOffsetAndSize(start, count);
    }
  else
    {
    offs_right |= data.M[3] >> 4;
    offs_right &= MAX_OFFSET_MASK;
    }
  }


/// Intialize from a "full" (not compressed) treeletized/clusterized tree.
/// @param[in] a_nodes        - array of nodes, may merge several BVH tree
/// @param[in] a_bvhOffsets   - array of positions (in a_nodes) where these trees start
/// @param[in] a_treelet_root - (indices of) the roots of the LAST LEVEL (smallest, 
///                             those which are not further clusterized) clusters.
///                             These are indices in the global array @c a_nodes
///                             (that may merge several BVH trees), i.e. do NOT add 'offset'.
/// @param[in] grStart        - (indices of) the starts of GROUPS of nodes.
///                             The bounding boxes for all nodes from one GROUP are 
///                             all calculated relative to the COMMON bounding box
///                             for the group (usually the union of the nodes B.B.)
///                             and in the array of compact storage is common B.B.
///                             is placed before the 1st node of that group.
///                             It is possible to use @c a_treelet_start for this argument
///                             then groups = treelets.
///                             These are indices in the global array @c a_nodes
///                             (that may merge several BVH trees), i.e. do NOT add 'offset'.
/// 
/// @return SUCCESS/FAILURE
int NonunifNodeStorage::Create(const std::vector<BVHNodeFat>& a_nodes, 
                               const std::vector<uint32_t>& a_bvhOffsets,
                               const std::vector<int>& a_treelet_root,
                               const std::vector<int>& grStart)
  {
  size_t i, j;
  // Check if all groups keep the same number of nodes and if it is 2^n-1
  uint32_t group_size = grStart.size() == 1 ? a_nodes.size() : grStart[1] - grStart[0];

  for (i = 1; i < grStart.size() - 1; i++)
    {
    if (grStart[i + 1] - grStart[i] != group_size)
      {
      assert(false);
      return -1;
      }
    }

  m_log2_group_size = -1;
  j = group_size + 1;
  for (i = 1; i < 1000; i++)
    {
    j = j >> 1;
    if (j == 1)
      {
      m_log2_group_size = i;
      break;
      }
    }
  if (m_log2_group_size <= 0 || group_size + 1 != (1 << m_log2_group_size))
    {
    assert(false);
    return -1;
    }
  //printf("%i nodes in treelet, log2(1 + ...) = %i\n", group_size, m_log2_group_size);




  // Size of the treelet's root B.B. data in the node array elements
  const int bb_block_size = TreeletBBoxTakesElements();

  const size_t n_nodes_all = a_nodes.size(), n_subroots = grStart.size();

  // Count the necessary number of 16 bytes blocks (to keep the data):
#ifdef HALFFLOAT
  const size_t n_16bytes = n_nodes_all + n_subroots;
#else
  const size_t n_16bytes = n_nodes_all + n_subroots * 2;
#endif

  data.resize(n_16bytes);

  // lookup from position of the 16bytes block in 'data' to the index of node (in a_nodes!) 
  // which starts at this position

  // lookup from the RELATIVE index of node (counted from tree's start in a_nodes!) 
  // to position of the 16bytes block in 'data' the compressed node data starts from
  std::vector<size_t> index2block;
  index2block.resize(n_nodes_all);

  m_bvhOffsets16B.resize(a_bvhOffsets.size());

  int curr_root_idx = -1; // index of the root of cluster current node belongs to
  int next_root_idx = -1; // index of the NEXT root of cluster
  int curr_root_number = -1; // The NUMBER of current root (=how many roots "before" it)

  m_bvhOffsets16B[0] = 0;
  for (j = 0; j < a_bvhOffsets.size(); j++)
    {
    const uint32_t offset = a_bvhOffsets[j];
    const size_t offset16B = m_bvhOffsets16B[j];

    const int n_nodes = j + 1 < a_bvhOffsets.size() ? a_bvhOffsets[j + 1] - offset : a_nodes.size() - offset;

    size_t curr_block = 0; // position in 'data' the current node must start at
    BBox3f curr_bbox; // box bounding left and right children of the current group
    LiteMath::uint4* nodeC; // compressed node data
    LiteMath::uint4* curr_root = NULL; // compressed node data for the starting node of the group

    LiteMath::float3 tr_bb_origin, tr_bb_ssize;// Actual parameters of the treelet's root BB (after its truncation for compact storage!)
    LiteMath::float3  tr_bbox_rev_ssize; // 1/(treelet's root bounding box scaled size)
    for (size_t i = 0; i < n_nodes; i++)
      {
      const BVHNodeFat& node = a_nodes[i + offset];

      const BBox3f cb1 = GetChildBoxLeft(node);
      const BBox3f cb2 = GetChildBoxRight(node);

      if (curr_root_number == -1 || i + offset == next_root_idx)
        {
        // next group begins here! 
        curr_root_number++;
        curr_root_idx = grStart[curr_root_number];
        next_root_idx = (curr_root_number + 1) < grStart.size() ? grStart[curr_root_number + 1] : a_nodes.size();

        // We shall first write two float which represent bounding box of the treelet's root 
        // so position for the node data shifts by two 16byte blocks

        curr_block += bb_block_size;

        curr_root = data.data() + curr_block + offset16B;

#if 1 
        // Common bounding box for this gropup of nodes 
        curr_bbox = CalcCommonBB(a_nodes, curr_root_idx, next_root_idx - 1);
#else
        // bounding box comprising cb1 and cb2
        curr_bbox.boxMin[0] = std::min(cb1.boxMin[0], cb2.boxMin[0]);
        curr_bbox.boxMin[1] = std::min(cb1.boxMin[1], cb2.boxMin[1]);
        curr_bbox.boxMin[2] = std::min(cb1.boxMin[2], cb2.boxMin[2]);
        curr_bbox.boxMax[0] = std::max(cb1.boxMax[0], cb2.boxMax[0]);
        curr_bbox.boxMax[1] = std::max(cb1.boxMax[1], cb2.boxMax[1]);
        curr_bbox.boxMax[2] = std::max(cb1.boxMax[2], cb2.boxMax[2]);
#endif


        assert(i + offset == curr_root_idx); // this node is just the root!

        SetRootBBox(curr_bbox, curr_root);

        // Retrieve what is ACTUALLY stored as the root BB parameters (as subjected to riundoff!)
        GetRootBBox(tr_bb_origin, tr_bb_ssize, curr_root);
        tr_bbox_rev_ssize = 1 / tr_bb_ssize;
        }
      assert(curr_bbox.boxMin.x >= -1e5 && curr_bbox.boxMin.x <= 1e5);
      assert(curr_bbox.boxMax.x >= -1e5 && curr_bbox.boxMax.x <= 1e5);
      assert(curr_bbox.boxMin.y >= -1e5 && curr_bbox.boxMin.y <= 1e5);
      assert(curr_bbox.boxMax.y >= -1e5 && curr_bbox.boxMax.y <= 1e5);
      assert(curr_bbox.boxMin.z >= -1e5 && curr_bbox.boxMin.z <= 1e5);
      assert(curr_bbox.boxMax.z >= -1e5 && curr_bbox.boxMax.z <= 1e5);

      // Conversion from the GLOBAL (with offset!) original node index to 
      // RELATIVE (from the start of the current tree) position of the node 
      // in the 'data' array
      index2block[i + offset] = curr_block; // "relative" offset, from this tree start!

      // compressed node data
      nodeC = data.data() + curr_block + offset16B;

      if (node.offs_left != 0xFFFFFFFF) // skip dummy node (added only for alignment)
        {
        // compressed child boxes
        SetChildBoxes(cb1, cb2, tr_bb_origin, tr_bbox_rev_ssize, nodeC);

#ifdef OSE_TESTS
        LiteMath::BBox3f cb1_, cb2_;
        GetChildBoxes(nodeC, cb1_, cb2_);
#endif
        }

      // After writing the node data which shifts position by one 16bytes block
      curr_block += 1;
      }

    if (j + 1 < m_bvhOffsets16B.size())
      m_bvhOffsets16B[j + 1] = curr_block + offset16B; // start of the NEXT tree will be here

    // Now that we know positions of all compressed nodes in the array, 
    // we can calculate what the original offs_left/offs_right and set them.
    for (size_t i = 0; i < n_nodes; i++)
      {
      const BVHNodeFat& node = a_nodes[i + offset];
      if (node.offs_left == 0xFFFFFFFF) // skip dummy node (added only for alignment)
        continue;

      uint offs_left, offs_right;
      if (!(node.offs_left & LEAF_BIT))
        offs_left = index2block[node.offs_left + offset]; /// TODO: Add the flag bits to the 25-bit representation!
      else
        offs_left = node.offs_left; /// TODO: Add the flag bits to the 25-bit representation!

      if (!(node.offs_right & LEAF_BIT))
        offs_right = index2block[node.offs_right + offset]; /// TODO: Add the flag bits to the 25-bit representation!
      else
        offs_right = node.offs_right; /// TODO: Add the flag bits to the 25-bit representation!

      nodeC = data.data() + index2block[i + offset] + offset16B;
      LiteMath::uint4 nodeC0 = *nodeC;

      SetOffsets(offs_left, offs_right, nodeC);

#ifdef OSE_TESTS
      uint offs_left1 = 0, offs_right1 = 0;
      GetOffsets(offs_left1, offs_right1, *nodeC);
      if (offs_left1 != offs_left || offs_right1 != offs_right)
        assert(false);
#endif
      }
    }

#ifdef OSE_TESTS
  for (j = 0; j < a_bvhOffsets.size(); j++)
    {
    const uint offset = a_bvhOffsets[j];
    const size_t offset16B = m_bvhOffsets16B[j];

    const int n_nodes = j + 1 < a_bvhOffsets.size() ? a_bvhOffsets[j + 1] - offset : a_nodes.size() - offset;

    //size_t curr_block = 0; // position in 'data' the current node must start at
    BBox3f curr_bbox; // box bounding left and right children of the current cluster's root
    LiteMath::uint4* nodeC; // compressed node data

    BBox3f cb1_, cb2_;

    for (size_t i = 0; i < n_nodes; i++)
      {
      const BVHNodeFat& node = a_nodes[i + offset];
      if (node.offs_left == 0xFFFFFFFF) // skip dummy node (added only for alignment)
        continue;

      const BBox3f cb1 = GetChildBoxLeft(node);
      const BBox3f cb2 = GetChildBoxRight(node);

      //curr_block = index2block[i + offset];

      nodeC = data.data() + index2block[i + offset] + offset16B;

      uint offs_left1 = 0, offs_right1 = 0;
      GetOffsets(offs_left1, offs_right1, *nodeC);

      if (!(node.offs_left & LEAF_BIT))
        {
        if (index2block[node.offs_left + offset] != offs_left1) 
          assert(false);
        }
      else
        {
        if (node.offs_left != offs_left1)
          assert(false);
        }
      if (!(node.offs_right & LEAF_BIT))
        {
        if (index2block[node.offs_right + offset] != offs_right1)
          assert(false);
        }
      else
        {
        if (node.offs_right != offs_right1)
          assert(false);
        }
      }
    }
#endif
  return 0;
  }

/// Appends the data, e.g. for a new geometry objected added to the scene
/// @param[in] nns --- BVH tree in this format to add to the current tree
///                    Unless the current object (='this') is empty,
///                    MUST BE OF THE SAME GROUP SIZE as 'this'!
/// @return SUCCESS/FAILURE
int NonunifNodeStorage::Append(const NonunifNodeStorage& nns)
  {
  const uint oldOffsetSize = m_bvhOffsets16B.size();
  const uint oldDataSize = data.size();

  if (oldDataSize > 0 && nns.m_log2_group_size != m_log2_group_size)
    {
    assert(false);
    return -1;
    }

  m_log2_group_size = nns.m_log2_group_size;

  data.insert(data.end(), nns.data.begin(), nns.data.end());
  m_bvhOffsets16B.insert(m_bvhOffsets16B.end(), nns.m_bvhOffsets16B.begin(), nns.m_bvhOffsets16B.end());

  if (oldDataSize > 0)
    {
    for (size_t i = oldOffsetSize; i < m_bvhOffsets16B.size(); i++)
      m_bvhOffsets16B[i] += oldDataSize;
    }
  return 0;
  }


/// Given the pointer to the compressed node data, calculates the bounding boxes of its left and right children
/// The input data can relate to the "internal" or "root" node.
/// It takes the pointer to the treelet root, shifts it left by 2*sizeof(float4) and
/// this is the pointer to the "pivot" bounding box stored as two float4.
/// Then it extracts from the current node 2x3 6bit numbers that represent
/// the left bounding box, RELATIVE to that of the root. Same for the right
/// bounding box. Scales then by the "pivot" box, and voila
/// 
/// @param[in]  data  - pointer to the (start of) compressed node data
/// @param[out] left  - bounding box of the left child of the above node
/// @param[out] right - bounding box of the left child of the above node
void NonunifNodeStorage::GetChildBoxes(const LiteMath::uint4* data, BBox3f& left, BBox3f& right) const 
  {
  // Pointer to the root data
  const LiteMath::uint4* treelet_root = GetTreeletRoot(data);

  // Bounding box information of the treelet root is represented by two float4
  // values, precceding 'treelet_root'. First is the origin, 2nd is scaled size

#ifdef HALFFLOAT
  LiteMath::float3 origin, ssize;
  const BoxScalerHalf* boxscaler = reinterpret_cast<const BoxScalerHalf*>(treelet_root) - 1;
  boxscaler->Get(origin, ssize);
#else
  // Pointer to the treelet's root bounding box origin
  const LiteMath::float4 &origin = *(reinterpret_cast<const LiteMath::float4*>(treelet_root) - 2);

  // Pointer to the bounding box scaled size
  const LiteMath::float4 &ssize = *(&origin + 1);
#endif

  // Coordinates of the bounding box of the left and of the right child, RELATIVE
  // to the B.B. of the treelet root (then its coordinates are in [0,1] and encoded
  // into 6 bits (i.e. in the range [0,63])
  LiteMath::uchar left6[6], right6[6];
  GetChildBoxes(data, left6, right6);

  left.boxMin[0] = origin[0] + left6[0] * ssize[0];
  left.boxMin[1] = origin[1] + left6[1] * ssize[1];
  left.boxMin[2] = origin[2] + left6[2] * ssize[2];

  left.boxMax[0] = origin[0] + left6[3] * ssize[0];
  left.boxMax[1] = origin[1] + left6[4] * ssize[1];
  left.boxMax[2] = origin[2] + left6[5] * ssize[2];


  right.boxMin[0] = origin[0] + right6[0] * ssize[0];
  right.boxMin[1] = origin[1] + right6[1] * ssize[1];
  right.boxMin[2] = origin[2] + right6[2] * ssize[2];

  right.boxMax[0] = origin[0] + right6[3] * ssize[0];
  right.boxMax[1] = origin[1] + right6[4] * ssize[1];
  right.boxMax[2] = origin[2] + right6[5] * ssize[2];
  }

/// Another variant, where the "pivot" bounding box (to scale the 6bit relative data)
/// is not found internally, but is provided explcitly as an argument.
/// @param[in] data  - node data
/// @param[in] bb_data - bounding box data (bb_data[0] is the origin and bb_data[1] is 1/size)
///                      for the root of the treelet/cluster this node belongs to.
///                      It can be taken as (const float4 *)(GetTreeletRoot(data)-2)
/// @param[out] left - bounding box of the left child
/// @param[out] right - bounding box of the right child
void NonunifNodeStorage::GetChildBoxes(const LiteMath::uint4* data, const LiteMath::float4 bb_data[2], BBox3f& left, BBox3f& right)
  {
  // Bounding box information of the treelet root is represented by two float4
  // values, precceding 'treelet_root'. First is the origin, 2nd is scaled size

  // Pointer to the treelet's root bounding box origin
  const LiteMath::float4 &origin = bb_data[0];

  // Pointer to the bounding box scaled size
  const LiteMath::float4 &ssize = bb_data[1];

  // Coordinates of the bounding box of the left and of the right child, RELATIVE
  // to the B.B. of the treelet root (then its coordinates are in [0,1] and encoded
  // into 6 bits (i.e. in the range [0,63])
  LiteMath::uchar left6[6], right6[6];
  GetChildBoxes(data, left6, right6);

  left.boxMin[0] = origin[0] + left6[0] * ssize[0];
  left.boxMin[1] = origin[1] + left6[1] * ssize[1];
  left.boxMin[2] = origin[2] + left6[2] * ssize[2];

  left.boxMax[0] = origin[0] + left6[3] * ssize[0];
  left.boxMax[1] = origin[1] + left6[4] * ssize[1];
  left.boxMax[2] = origin[2] + left6[5] * ssize[2];


  right.boxMin[0] = origin[0] + right6[0] * ssize[0];
  right.boxMin[1] = origin[1] + right6[1] * ssize[1];
  right.boxMin[2] = origin[2] + right6[2] * ssize[2];

  right.boxMax[0] = origin[0] + right6[3] * ssize[0];
  right.boxMax[1] = origin[1] + right6[4] * ssize[1];
  right.boxMax[2] = origin[2] + right6[5] * ssize[2];
  }

/// Again the "pivot" bounding box (to scale the 6bit relative data)
/// is not found internally, but is provided explcitly as an argument, but as float3.
/// @param[in] data  - node data
/// @param[in] bb_origin - min point of the bounding box of the treelet's root
/// @param[in] bb_ssize  - size/63 of the bounding box of the treelet's root
/// @param[out] left - bounding box of the left child
/// @param[out] right - bounding box of the right child
void NonunifNodeStorage::GetChildBoxes(const LiteMath::uint4* data, const LiteMath::float3& bb_origin, const LiteMath::float3& bb_ssize, BBox3f& left, BBox3f& right)
  {
  // Coordinates of the bounding box of the left and of the right child, RELATIVE
  // to the B.B. of the treelet root (then its coordinates are in [0,1] and encoded
  // into 6 bits (i.e. in the range [0,63])
  LiteMath::uchar left6[6], right6[6];
  GetChildBoxes(data, left6, right6);

  left.boxMin[0] = bb_origin[0] + left6[0] * bb_ssize[0];
  left.boxMin[1] = bb_origin[1] + left6[1] * bb_ssize[1];
  left.boxMin[2] = bb_origin[2] + left6[2] * bb_ssize[2];

  left.boxMax[0] = bb_origin[0] + left6[3] * bb_ssize[0];
  left.boxMax[1] = bb_origin[1] + left6[4] * bb_ssize[1];
  left.boxMax[2] = bb_origin[2] + left6[5] * bb_ssize[2];


  right.boxMin[0] = bb_origin[0] + right6[0] * bb_ssize[0];
  right.boxMin[1] = bb_origin[1] + right6[1] * bb_ssize[1];
  right.boxMin[2] = bb_origin[2] + right6[2] * bb_ssize[2];

  right.boxMax[0] = bb_origin[0] + right6[3] * bb_ssize[0];
  right.boxMax[1] = bb_origin[1] + right6[4] * bb_ssize[1];
  right.boxMax[2] = bb_origin[2] + right6[5] * bb_ssize[2];
  }