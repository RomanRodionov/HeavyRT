#include "BVH2Stackless.h"

void BVH2Stackless::IntersectAllPrimitivesInLeaf(const float3 ray_pos, const float3 ray_dir,
                                                float tNear, uint32_t instId, uint32_t geomId,
                                                uint32_t a_start, uint32_t a_count,
                                                CRT_Hit *pHit)
{
  const uint2 a_geomOffsets = m_geomOffsets[geomId];

  for (uint32_t triId = a_start; triId < a_start + a_count; triId++)
  {
    const uint32_t A = m_indices[a_geomOffsets.x + triId*3 + 0];
    const uint32_t B = m_indices[a_geomOffsets.x + triId*3 + 1];
    const uint32_t C = m_indices[a_geomOffsets.x + triId*3 + 2];

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

    if (v >= -1e-6f && u >= -1e-6f && (u + v <= 1.0f + 1e-6f) && t > tNear && t < pHit->t) 
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

CRT_Hit BVH2Stackless::RayQuery_NearestHit(float4 posAndNear, float4 dirAndFar)
{
  bool stopOnFirstHit = (dirAndFar.w <= 0.0f);
  if(stopOnFirstHit)
    dirAndFar.w *= -1.0f;

  #ifdef ENABLE_METRICS
  ResetVarLC();
  m_stats.raysNumber++;
  #endif

  // (2) process all intersected BLAS
  //
  CRT_Hit hit;
  hit.t      = dirAndFar.w;
  hit.primId = uint32_t(-1);
  hit.instId = uint32_t(-1);
  hit.geomId = uint32_t(-1);
  
  ////////////////////////////////////////////////////////////////////////////// instancing variables are a bit more than common
  uint32_t bvhOffset = m_tlasOffset;
  uint32_t instId    = uint32_t(-1);
  float3 ray_pos     = to_float3(posAndNear);
  float3 ray_dir     = to_float3(dirAndFar);
  float3 inv_dir     = SafeInverse(ray_dir);
  ////////////////////////////////////////////////////////////////////////////// instancing variables are a bit more than common
  
  uint32_t nodeIdx    = 0;
  uint32_t nodeIdxTop = 0;
  
  do
  {
    bool isLeaf         = false;
    bool intersects     = false;
    uint32_t leftOffset = 0;

    //do                                                           // NOTE: uncomment this to get while-while traversal!
    {
      const BVHNode currNode = m_allNodes[bvhOffset + nodeIdx];
      const float2 boxHit    = RayBoxIntersection2(ray_pos, inv_dir, currNode.boxMin, currNode.boxMax);
      
      #ifdef ENABLE_METRICS
      m_stats.NC  += 1;
      m_stats.BLB += 1 * sizeof(BVHNode);
      #endif
     
      intersects = (boxHit.x <= boxHit.y) && (boxHit.y > posAndNear.w) && (boxHit.x < hit.t); // (tmin <= tmax) && (tmax > 0.f) && (tmin < curr_t)
      isLeaf     = ((currNode.leftOffset & LEAF_BIT) != 0);
     
      leftOffset = currNode.leftOffset;
      nodeIdx    = (isLeaf || !intersects) ? currNode.escapeIndex : leftOffset;
      
    } // while (!isLeaf && intersects && nodeIdx < 0xFFFFFFFE);    // NOTE: uncomment this to get while-while traversal!

    if(intersects && isLeaf && bvhOffset != m_tlasOffset) //leaf
    {
      const uint32_t start  = EXTRACT_START(leftOffset);
      const uint32_t count  = EXTRACT_COUNT(leftOffset);
      const uint32_t geomId = m_geomIdByInstId[instId];
      
      IntersectAllPrimitivesInLeaf(ray_pos, ray_dir, posAndNear.w, instId, geomId, start, count, &hit); //

      #ifdef ENABLE_METRICS
      m_stats.LC++;
      m_stats.LC2++;
      m_stats.TC+=count;
      m_stats.BLB+=count*(3*sizeof(uint32_t) + 3*sizeof(float4));
      #endif
    }

    if(isLeaf && intersects && bvhOffset == m_tlasOffset) 
    {
      instId     = EXTRACT_START(leftOffset);
      bvhOffset  = m_bvhOffsets[m_geomIdByInstId[instId]];
      nodeIdxTop = nodeIdx;
      nodeIdx    = 0;

      ray_pos = matmul4x3(m_instMatricesInv[instId], to_float3(posAndNear));
      ray_dir = matmul3x3(m_instMatricesInv[instId], to_float3(dirAndFar));
      inv_dir = SafeInverse(ray_dir);
    }
  
    if (nodeIdx >= 0xFFFFFFFE && bvhOffset != m_tlasOffset)
    {
      ray_pos   = to_float3(posAndNear);
      ray_dir   = to_float3(dirAndFar);
      inv_dir   = SafeInverse(ray_dir);
      bvhOffset = m_tlasOffset;
      nodeIdx   = nodeIdxTop;
    }

  } while (nodeIdx < 0xFFFFFFFE && !(stopOnFirstHit && hit.primId != uint32_t(-1)));
  
  #ifdef REMAP_PRIM_ID
  if(hit.geomId < uint32_t(-1)) 
  {
    const uint2 geomOffsets = m_geomOffsets[hit.geomId];
    hit.primId = m_primIndices[geomOffsets.x/3 + hit.primId];
  }
  #endif

  return hit;
}

bool BVH2Stackless::RayQuery_AnyHit(float4 posAndNear, float4 dirAndFar)
{
  dirAndFar.w *= -1.0f;
  CRT_Hit hit = RayQuery_NearestHit(posAndNear, dirAndFar);
  return (hit.geomId != uint32_t(-1));
}
