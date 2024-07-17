#include "ShtRT64.h"

static inline float2 RayBoxIntersection(float3 rayOrigin, float3 rayDirInv, float3 boxMin, float3 boxMax)
{
  const float lo = rayDirInv.x * (boxMin.x - rayOrigin.x);
  const float hi = rayDirInv.x * (boxMax.x - rayOrigin.x);
  const float lo1 = rayDirInv.y * (boxMin.y - rayOrigin.y);
  const float hi1 = rayDirInv.y * (boxMax.y - rayOrigin.y);
  const float lo2 = rayDirInv.z * (boxMin.z - rayOrigin.z);
  const float hi2 = rayDirInv.z * (boxMax.z - rayOrigin.z);

  const float tmin = std::max(std::min(lo, hi), std::min(lo1, hi1));
  const float tmax = std::min(std::max(lo, hi), std::max(lo1, hi1));

  return float2(std::max(tmin, std::min(lo2, hi2)),
                std::min(tmax, std::max(lo2, hi2)));

}


bool ShtRT64::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
#ifdef ENABLE_METRICS
  m_stats.raysNumber++;
#endif

  dirAndFar.w *= -1.0f;
  CRT_Hit hit = RayQuery_NearestHit(posAndNear, dirAndFar);
  return (hit.geomId != uint32_t(-1));
}

uint32_t ShtRT64::CountTrailingZeros(uint64_t bitTrail)
{
  uint32_t bits = 0;
  uint64_t x = bitTrail;

  if(x != 0)
  {
    ///* assuming `x` has 64 bits: lets count the low order 0 bits in batches */
    if((x & 0xFFFFFFFF) == 0)
    {
     //std::cout << "HIT" << std::endl;
     bits += 32;
     x >>= 32;
    }

    /* mask the 16 low order bits, add 16 and shift them out if they are all 0 */

    if((x & 0x0000FFFF) == 0)
    {
      bits += 16;
      x >>= 16;
    }
    /* mask the 8 low order bits, add 8 and shift them out if they are all 0 */
    if((x & 0x000000FF) == 0)
    {
      bits += 8;
      x >>= 8;
    }
    /* mask the 4 low order bits, add 4 and shift them out if they are all 0 */
    if((x & 0x0000000F) == 0)
    {
      bits += 4;
      x >>= 4;
    }
    /* mask the 2 low order bits, add 2 and shift them out if they are all 0 */
    if((x & 0x00000003) == 0)
    {
      bits += 2;
      x >>= 2;
    }
    /* mask the low order bit and add 1 if it is 0 */
    bits += (x & 1) ^ 1;
  }
  return bits;
}

void ShtRT64::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir, float tNear, uint32_t instId,
                                         uint32_t geomId, uint32_t a_start, uint32_t a_count, CRT_Hit *pHit)
{
  const uint2 a_geomOffsets = m_geomOffsets[geomId];

  for(uint32_t triId = a_start; triId < a_start + a_count; triId++)
  {
    const uint32_t A = m_indices[a_geomOffsets.x + triId * 3 + 0];
    const uint32_t B = m_indices[a_geomOffsets.x + triId * 3 + 1];
    const uint32_t C = m_indices[a_geomOffsets.x + triId * 3 + 2];

    const float3 A_pos = to_float3(m_vertPos[a_geomOffsets.y + A]);
    const float3 B_pos = to_float3(m_vertPos[a_geomOffsets.y + B]);
    const float3 C_pos = to_float3(m_vertPos[a_geomOffsets.y + C]);

    const float3 edge1 = B_pos - A_pos;
    const float3 edge2 = C_pos - A_pos;
    const float3 pvec = cross(ray_dir, edge2);
    const float3 tvec = ray_pos - A_pos;
    const float3 qvec = cross(tvec, edge1);

    const float invDet = 1.0f / dot(edge1, pvec);
    const float v = dot(tvec, pvec) * invDet;
    const float u = dot(qvec, ray_dir) * invDet;
    const float t = dot(edge2, qvec) * invDet;

    if(v >= -1e-6f && u >= -1e-6f && (u + v <= 1.0f + 1e-6f) && t > tNear && t < pHit->t) // if (v > -1e-6f && u > -1e-6f && (u + v < 1.0f+1e-6f) && t > tMin && t < hit.t)
    {
      pHit->t = t;
      pHit->primId = triId;
      pHit->instId = instId;
      pHit->geomId = geomId;
      pHit->coords[0] = u;
      pHit->coords[1] = v;
    }
  }
}

CRT_Hit ShtRT64::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  const bool stopOnFirstHit = (dirAndFar.w <= 0.0f);
  if(stopOnFirstHit)
    dirAndFar.w *= -1.0f;
    
  #ifdef ENABLE_METRICS
  ResetVarLC();
  m_stats.raysNumber++;
  #endif
  
  CRT_Hit hit;
  hit.t = dirAndFar.w;
  hit.primId = uint32_t(-1);
  hit.instId = uint32_t(-1);
  hit.geomId = uint32_t(-1);

  const float3 rayDirInv = SafeInverse(to_float3(dirAndFar));
  uint32_t nodeIdx = 0;
  do
  {
    uint32_t travFlags  = 0;
    uint32_t leftOffset = 0;
    do
    {
      const BVHNode currNode = m_nodesTLAS[nodeIdx];
      const float2 boxHit    = RayBoxIntersection2(to_float3(posAndNear), rayDirInv, currNode.boxMin, currNode.boxMax);
      const bool intersects  = (boxHit.x <= boxHit.y) && (boxHit.y > posAndNear.w) && (boxHit.x < hit.t); // (tmin <= tmax) && (tmax > 0.f) && (tmin < curr_t)
      
      #ifdef ENABLE_METRICS
      m_stats.NC  += 1;
      m_stats.BLB += 1 * sizeof(BVHNode);
      #endif

      travFlags  = (currNode.leftOffset & LEAF_BIT) | uint32_t(intersects); // travFlags  = (((currNode.leftOffset & LEAF_BIT) == 0) ? 0 : LEAF_BIT) | (intersects ? 1 : 0);
      leftOffset = currNode.leftOffset;
      nodeIdx    = isLeafOrNotIntersect(travFlags) ? currNode.escapeIndex : leftOffset;

    } while (notLeafAndIntersect(travFlags) && nodeIdx != 0 && nodeIdx < 0xFFFFFFFE); 
     
    if(isLeafAndIntersect(travFlags)) 
    {
      const uint32_t instId = EXTRACT_START(leftOffset);
      const uint32_t geomId = m_geomIdByInstId[instId];
  
      // transform ray with matrix to local space
      //
      const float3 ray_pos = matmul4x3(m_instMatricesInv[instId], to_float3(posAndNear));
      const float3 ray_dir = matmul3x3(m_instMatricesInv[instId], to_float3(dirAndFar)); // DON'T NORMALIZE IT !!!! When we transform to local space of node, ray_dir must be unnormalized!!!
      Traverse(ray_pos, ray_dir, posAndNear.w, instId, geomId, &hit);
    }

  } while (nodeIdx < 0xFFFFFFFE && !(stopOnFirstHit && hit.primId != uint32_t(-1)));

  #ifdef REMAP_PRIM_ID
  if(hit.geomId < uint32_t(-1))
  {
    const uint2 geomOffsets = m_geomOffsets[hit.geomId];
    hit.primId = m_primIndices[geomOffsets.x / 3 + hit.primId];
  }
  #endif
  
  return hit;
}

void ShtRT64::Traverse(const float3 rayPos, const float3 ray_dir, float tMin, uint32_t instId, uint32_t geomId, CRT_Hit *hit)
{
  const uint32_t bvhOffset = m_bvhOffsets[geomId];

  uint32_t currentNode = 0;
  uint32_t oldLeftNodeOffset = 0;
  uint32_t leftNodeOffset = 0;
  uint64_t nodeKey = 1; //key for root; for node with key=k its left_key=2 * k and right_key=2 * k + 1
  uint64_t bitTrail = 0;
  uint32_t mrpnAddr = 0xFFFFFFFF; // most recently postponed node
  const float3 rayDirInv = SafeInverse(ray_dir);
  #ifdef LINEAR_HASHING
  uint32_t hOffset = hOffsets[geomId];
  uint32_t dOffset = dOffsets[geomId];
  uint32_t dSize = dispSizes[geomId];
  uint32_t hSize = hSizes[geomId];
  #endif

  #ifdef LINEAR_HASHING
  #else
  HashTable *hashTable = hashTables[geomId];
  #endif

  #ifdef ENABLE_METRICS
  uint32_t oldCurrentNode = currentNode;
  #endif

  while(oldLeftNodeOffset != 0xFFFFFFFF)
  { 
    uint32_t escapeIndex;
    { 
      const BVHNode oldNode = m_allNodes[bvhOffset + oldLeftNodeOffset];
      currentNode    = oldLeftNodeOffset;
      leftNodeOffset = oldNode.leftOffset;
      escapeIndex        = oldNode.escapeIndex;
    }

    //forward Traverse until find leaf
    bool searchingForLeaf = (leftNodeOffset & LEAF_BIT) == 0 && (leftNodeOffset != EMPTY_NODE);
    while(searchingForLeaf) //root node index
    {
      #ifdef ENABLE_METRICS
      m_stats.NC  += 2;
      m_stats.BLB += 2 * sizeof(BVHNode);
      #endif
      const BVHNode leftNode  = m_allNodes[bvhOffset + leftNodeOffset + 0];
      const BVHNode rightNode = m_allNodes[bvhOffset + leftNodeOffset + 1];
      const uint2 nextOffsets = uint2(leftNode.leftOffset, rightNode.leftOffset);

      //TODO add uncle, grand uncle
      if(mrpnAddr == 0xFFFFFFFF && (bitTrail & 0x2) != 0)
        mrpnAddr = escapeIndex;

      const float2 t_left  = RayBoxIntersection(rayPos, rayDirInv, leftNode.boxMin, leftNode.boxMax);
      const float2 t_right = RayBoxIntersection(rayPos, rayDirInv, rightNode.boxMin, rightNode.boxMax);

      const bool isHitLeftChild  = (t_left.x  <= t_left.y)  && (t_left.y >= tMin)  && (t_left.x <= hit->t);
      const bool isHitRightChild = (t_right.x <= t_right.y) && (t_right.y >= tMin) && (t_right.x <= hit->t); //bad t_right for 88303
      const bool closetIsLeft    = (t_left.x  <= t_right.x) || !isHitRightChild;
      
      //identify first hit
      if(isHitLeftChild || isHitRightChild)
      {
        //#ifdef ENABLE_METRICS
        //if (isHitLeftChild && isHitRightChild)
        //  m_stats.SOC++;
        //#endif
        currentNode = (isHitLeftChild && closetIsLeft) ? leftNodeOffset : leftNodeOffset + 1;
        nodeKey     = nodeKey << 1;
        bitTrail    = bitTrail << 1;

        if(isHitLeftChild && isHitRightChild)
        {
          bitTrail = bitTrail ^ 1;
          mrpnAddr = closetIsLeft ? leftNodeOffset + 1 : leftNodeOffset;
        }
        if(currentNode == leftNodeOffset + 1)
          nodeKey = nodeKey ^ 1;
        
        leftNodeOffset = (isHitLeftChild && closetIsLeft) ? nextOffsets.x : nextOffsets.y;
      }
     
      searchingForLeaf = (leftNodeOffset & LEAF_BIT) == 0 && (leftNodeOffset != EMPTY_NODE) && (isHitLeftChild || isHitRightChild);
    }

    // leaf node found
    if((leftNodeOffset & LEAF_BIT) != 0 && leftNodeOffset != 0xFFFFFFFF)
    {
      // intersect
      auto start = EXTRACT_START(leftNodeOffset);
      auto count = EXTRACT_COUNT(leftNodeOffset);
      IntersectAllPrimitivesInLeaf(rayPos, ray_dir, tMin, instId, geomId, start, count, hit);
      if(mrpnAddr == 0xFFFFFFFF && (bitTrail & 0x2) != 0)
        mrpnAddr = escapeIndex;

      #ifdef ENABLE_METRICS
      m_stats.LC++;
      m_stats.LC2++;
      m_stats.TC += count;
      m_stats.BLB += count * (3 * sizeof(uint32_t) + 3 * sizeof(float4));
      #endif
    }

    if(bitTrail == 0)
    {
      // No backtracking
      oldLeftNodeOffset = 0xFFFFFFFF;
    }
    else
    {
      // Backtracking
      uint32_t numLevels = CountTrailingZeros(bitTrail);
      bitTrail = (bitTrail >> numLevels) ^ 1;
      nodeKey = (nodeKey >> numLevels) ^ 1;

      if(mrpnAddr != 0xFFFFFFFF)
      {
        oldLeftNodeOffset = mrpnAddr;
        mrpnAddr = 0xFFFFFFFF;
      }
      else
      {
        //TODO always 0 for NOW
        #ifdef LINEAR_HASHING
        oldLeftNodeOffset = hashTable[hOffset + (nodeKey + displacementTables[dOffset + (nodeKey & (dSize - 1))]) % hSize];
        #else
        oldLeftNodeOffset = hashTable->get(nodeKey);
        #endif
      }

      #ifdef ENABLE_METRICS
      m_stats.SOC++;
      #endif
    }
  }
}

