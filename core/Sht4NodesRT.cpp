#include "Sht4NodesRT.h"

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

bool Sht4NodesRT::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
#ifdef ENABLE_METRICS
    m_stats.raysNumber++;
#endif

    // (2) If any hit is found, immediately return true.
    std::cout << "[Sht4NodesRT::RayQuery_AnyHit]: "
              << "not implemeted!" << std::endl;

    return false;
}

uint32_t Sht4NodesRT::CountTrailingZeros(uint64_t bitTrail)
{
    uint64_t bits = 0, x = bitTrail;

    if(x != 0)
    {
        ///* assuming `x` has 64 bits: lets count the low order 0 bits in batches */
        if((x & 0xFFFFFFFF) == 0)
        {
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

void Sht4NodesRT::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir, float tNear, uint32_t instId,
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


CRT_Hit Sht4NodesRT::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  #ifdef ENABLE_METRICS
  m_stats.raysNumber++;
  ResetVarLC();
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

  } while (nodeIdx < 0xFFFFFFFE);

  #ifdef REMAP_PRIM_ID
  if(hit.geomId < uint32_t(-1))
  {
    const uint2 geomOffsets = m_geomOffsets[hit.geomId];
    hit.primId = m_primIndices[geomOffsets.x / 3 + hit.primId];
  }
  #endif
  return hit;
}

uint32_t num_c = 0;

void
Sht4NodesRT::Traverse(const float3 rayPos, const float3 ray_dir, float tMin, uint32_t instId, uint32_t geomId,
                      CRT_Hit *hit)
{
    const uint32_t bvhOffset = m_bvhOffsets[geomId];
    uint32_t leftNodeOffset = 0;
    #ifdef ENABLE_METRICS
    auto leftNodeOffsetOld2 = leftNodeOffset;
    #endif
    int top = 0;
    uint32_t stack[6];
    uint64_t nodeKey = 1; //key for root; for node with key=k its left_key=2 * k and right_key=2 * k + 1
    uint64_t postponedNodes = 0;
    uint64_t leftOnLevel = 0;
    uint64_t numLevels = 0;
    uint64_t zeroTwoBits = 0xFFFFFFFFFFFFFFFC;

    const float3 rayDirInv = SafeInverse(ray_dir);
    uint32_t hOffset = hOffsets[geomId];
    uint32_t dOffset = dOffsets[geomId];
    uint32_t dSize = dispSizes[geomId] - 1;
    uint32_t hSize = hSizes[geomId];

    while(true)
    {
        //forward Traverse until leaf
        bool searchingForLeaf = (leftNodeOffset & LEAF_BIT) == 0 &&
                                (leftNodeOffset != 0xFFFFFFFF);
        while (searchingForLeaf)
        {
#ifdef ENABLE_METRICS
            m_stats.NC += 4;
            m_stats.BLB += 4 * sizeof(BVHNode);
            if(leftNodeOffsetOld2 != leftNodeOffset)
            {
                for(int i = 0; i < TREELET_ARR_SIZE; i++)
                {
                    if(std::abs(double(leftNodeOffset) - double(leftNodeOffsetOld2)) * double(sizeof(BVHNode)) >=
                       (double) treelet_sizes[i])
                        m_stats.LJC[i]++;
                    const uint32_t oldCacheLineId =
                            uint32_t(leftNodeOffsetOld2 * sizeof(BVHNode)) / uint32_t(treelet_sizes[i]);
                    const uint32_t newCacheLineId =
                            uint32_t(leftNodeOffset * sizeof(BVHNode)) / uint32_t(treelet_sizes[i]);
                    if(oldCacheLineId != newCacheLineId)
                    {
                        m_stats.CMC[i]++;
                        m_stats.WSS[i].insert(newCacheLineId);
                    }
                }
                leftNodeOffsetOld2 = leftNodeOffset;
            }
#endif
            const BVHNode node0 = m_allNodes[bvhOffset + leftNodeOffset + 0]; 
            const BVHNode node1 = m_allNodes[bvhOffset + leftNodeOffset + 1]; 
            const BVHNode node2 = m_allNodes[bvhOffset + leftNodeOffset + 2]; 
            const BVHNode node3 = m_allNodes[bvhOffset + leftNodeOffset + 3]; 

            uint4 leftOffsets = uint4(node0.leftOffset, node1.leftOffset, node2.leftOffset, node3.leftOffset);

            const float2 tm0 = RayBoxIntersection2(rayPos, rayDirInv, node0.boxMin, node0.boxMax);
            const float2 tm1 = RayBoxIntersection2(rayPos, rayDirInv, node1.boxMin, node1.boxMax);
            const float2 tm2 = RayBoxIntersection2(rayPos, rayDirInv, node2.boxMin, node2.boxMax);
            const float2 tm3 = RayBoxIntersection2(rayPos, rayDirInv, node3.boxMin, node3.boxMax);

            const bool hitChild0 = (tm0.x <= tm0.y) && (tm0.y >= tMin) && (tm0.x <= hit->t);
            const bool hitChild1 = (tm1.x <= tm1.y) && (tm1.y >= tMin) && (tm1.x <= hit->t);
            const bool hitChild2 = (tm2.x <= tm2.y) && (tm2.y >= tMin) && (tm2.x <= hit->t);
            const bool hitChild3 = (tm3.x <= tm3.y) && (tm3.y >= tMin) && (tm3.x <= hit->t);
            const int toPostpone = hitChild0 + hitChild1 + hitChild2 + hitChild3 - 1;
            uint4 oldPositions = uint4(0,1,2,3);
            float4 hitMinD = float4(hitChild0 ? tm0.x : MAXFLOAT,
                                    hitChild1 ? tm1.x : MAXFLOAT,
                                    hitChild2 ? tm2.x : MAXFLOAT,
                                    hitChild3 ? tm3.x : MAXFLOAT);

            // sort tHit and oldPositions
            //
            const bool lessXY = (hitMinD.y < hitMinD.x);
            const bool lessWZ = (hitMinD.w < hitMinD.z);
            {
                const int   w_oldPositionsX = lessXY ? oldPositions.y : oldPositions.x;
                const int   w_oldPositionsY = lessXY ? oldPositions.x : oldPositions.y;
                
                const int   w_leftOffsetX = lessXY ? leftOffsets.y : leftOffsets.x;
                const int   w_leftOffsetY = lessXY ? leftOffsets.x : leftOffsets.y;
                
                const float w_hitTimesX = lessXY ? hitMinD.y  : hitMinD.x;
                const float w_hitTimesY = lessXY ? hitMinD.x  : hitMinD.y;

                const int   w_leftOffsetZ = lessWZ ? leftOffsets.w : leftOffsets.z;
                const int   w_leftOffsetW = lessWZ ? leftOffsets.z : leftOffsets.w;
                
                const int   w_oldPositionsZ = lessWZ ? oldPositions.w : oldPositions.z;
                const int   w_oldPositionsW = lessWZ ? oldPositions.z : oldPositions.w;

                const float w_hitTimesZ = lessWZ ? hitMinD.w  : hitMinD.z;
                const float w_hitTimesW = lessWZ ? hitMinD.z  : hitMinD.w;

                leftOffsets.x = w_leftOffsetX;
                leftOffsets.y = w_leftOffsetY;
                
                oldPositions.x = w_oldPositionsX;
                oldPositions.y = w_oldPositionsY;
                
                hitMinD.x  = w_hitTimesX;
                hitMinD.y  = w_hitTimesY;

                leftOffsets.z = w_leftOffsetZ;
                leftOffsets.w = w_leftOffsetW;

                oldPositions.z = w_oldPositionsZ;
                oldPositions.w = w_oldPositionsW;
                hitMinD.z  = w_hitTimesZ;
                hitMinD.w  = w_hitTimesW;
            }

            const bool lessZX = (hitMinD.z < hitMinD.x);
            const bool lessWY = (hitMinD.w < hitMinD.y);
            {
                const int   w_leftOffsetX = lessZX ? leftOffsets.z : leftOffsets.x;
                const int   w_leftOffsetZ = lessZX ? leftOffsets.x : leftOffsets.z;

                const int   w_oldPositionsX = lessZX ? oldPositions.z : oldPositions.x;
                const int   w_oldPositionsZ = lessZX ? oldPositions.x : oldPositions.z;
                const float w_hitTimesX = lessZX ? hitMinD.z  : hitMinD.x;
                const float w_hitTimesZ = lessZX ? hitMinD.x  : hitMinD.z;

                const int   w_leftOffsetY = lessWY ? leftOffsets.w : leftOffsets.y;
                const int   w_leftOffsetW = lessWY ? leftOffsets.y : leftOffsets.w;
                const int   w_oldPositionsY = lessWY ? oldPositions.w : oldPositions.y;
                const int   w_oldPositionsW = lessWY ? oldPositions.y : oldPositions.w;
                const float w_hitTimesY = lessWY ? hitMinD.w  : hitMinD.y;
                const float w_hitTimesW = lessWY ? hitMinD.y  : hitMinD.w;

                leftOffsets.x = w_leftOffsetX;
                leftOffsets.z= w_leftOffsetZ;
                oldPositions.x = w_oldPositionsX;
                oldPositions.z = w_oldPositionsZ;
                hitMinD.x  = w_hitTimesX;
                hitMinD.z  = w_hitTimesZ;

                leftOffsets.y = w_leftOffsetY;
                leftOffsets.w= w_leftOffsetW;
                oldPositions.y = w_oldPositionsY;
                oldPositions.w = w_oldPositionsW;
                hitMinD.y  = w_hitTimesY;
                hitMinD.w  = w_hitTimesW;
            }

            const bool lessZY = (hitMinD.z < hitMinD.y);
            {
                const int   w_oldPositionsY = lessZY ? oldPositions.z : oldPositions.y;
                const int   w_oldPositionsZ = lessZY ? oldPositions.y : oldPositions.z;
                const int   w_leftOffsetsY = lessZY ? leftOffsets.z : leftOffsets.y;
                const int   w_leftOffsetsZ = lessZY ? leftOffsets.y : leftOffsets.z;
                const float w_hitTimesY = lessZY ? hitMinD.z  : hitMinD.y;
                const float w_hitTimesZ = lessZY ? hitMinD.y  : hitMinD.z;

                leftOffsets.y = w_leftOffsetsY;
                leftOffsets.z = w_leftOffsetsZ;
                oldPositions.y = w_oldPositionsY;
                oldPositions.z = w_oldPositionsZ;
                hitMinD.y  = w_hitTimesY;
                hitMinD.z  = w_hitTimesZ;
            }

            if (toPostpone > 0) {
              postponedNodes <<= toPostpone * 2;
              leftOnLevel <<= 2;
              leftOnLevel |= toPostpone;
              numLevels <<= 1;
              numLevels |= 1;
            }
            if(hitMinD.w < MAXFLOAT)
            {
                postponedNodes |= oldPositions[3] << 4;
            }
            if(hitMinD.z < MAXFLOAT)
            {
                postponedNodes |= oldPositions[2] << 2;
            }
            if(hitMinD.y < MAXFLOAT)
            {
                postponedNodes |= oldPositions[1];
            }
            if(hitMinD.x < MAXFLOAT)
            {
                nodeKey = (nodeKey << 2) + oldPositions[0];
                if (hitMinD.y >= MAXFLOAT) {
                  numLevels <<= 1;
                }
            }
            else
            {
                break;
            }
            leftNodeOffset = leftOffsets.x;
            searchingForLeaf = (leftNodeOffset & LEAF_BIT) == 0 &&
                                (leftNodeOffset != 0xFFFFFFFF);
        }
        // leaf node found
        if((leftNodeOffset & LEAF_BIT) != 0 && leftNodeOffset < 0xFFFFFFFD)
        {
            // intersect
            const uint32_t start = EXTRACT_START(leftNodeOffset);
            const uint32_t count = EXTRACT_COUNT(leftNodeOffset);
            IntersectAllPrimitivesInLeaf(rayPos, ray_dir, tMin, instId, geomId, start, count,
                                         hit);
            #ifdef ENABLE_METRICS
            m_stats.LC++;
            m_stats.LC2++;
            m_stats.TC += count;
            m_stats.BLB += count * (3 * sizeof(uint32_t) + 3 * sizeof(float4));
            #endif
        }

        if (numLevels == 0) break;
        
        uint32_t numLev = CountTrailingZeros(numLevels);
        numLevels >>= numLev; 
        uint32_t pos = postponedNodes & 0x3;
        postponedNodes >>= 2;
        leftOnLevel--;
        if ((leftOnLevel & 0x3) == 0) {
            leftOnLevel >>= 2;
            numLevels ^= 1;
        }
        nodeKey = ((nodeKey >> (2 * numLev)) & zeroTwoBits) + pos;
        leftNodeOffset = hashTable[hOffset + (nodeKey + displacementTables[dOffset + (nodeKey & dSize)]) % hSize];
    }
}


void Sht4NodesRT::UpdateInstance(uint32_t a_instanceId, const float4x4 &a_matrix)
{
  std::cout << "[Sht4NodesRT::UpdateInstance]: "
            << "not implemeted!" << std::endl;
}


