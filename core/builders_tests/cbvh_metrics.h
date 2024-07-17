#pragma once
#include "cbvh_core.h"
#include "cbvh_internal.h"

using namespace cbvh_internal;

namespace cbvh_metrics{

  /**
  \brief Calculate SAH metric for BVH tree

  \param  in_bvhTree - input bvh tree
  \param  in_verts   - input triangle vertices
  \param  in_innerC  - input node cost
  \param  in_leafC   - input leaf cost

  \return metric's value
  */
  float MetricSAH(const cbvh_internal::BVHTree &in_bvhTree, const float4* in_verts, unsigned int in_verts_count,
                  const float in_innerC = 1.0, const float in_leafC = 1.0);

  /**
  \brief Calculate SAH metric for BVH tree using triangle-boxes for calculation of overlapping

  \param  in_bvhTree - input bvh tree
  \param  in_verts   - input triangle vertices
  \param  in_innerC  - input node cost
  \param  in_leafC   - input leaf cost

  \return metric's value
  */
  float MetricEPOTriBF(const cbvh_internal::BVHTree &in_bvhTree, const float4* in_verts, unsigned int in_verts_count,
                    const float in_innerC = 1.0, const float in_leafC = 1.0);

  /**
  \brief Calculate SAH metric for BVH tree using triangle-boxes for calculation of overlapping

  \param  in_bvhTree - input bvh tree
  \param  in_verts   - input triangle vertices
  \param  in_innerC  - input node cost
  \param  in_leafC   - input leaf cost

  \return metric's value
  */
  float MetricEPOTri(const cbvh_internal::BVHTree &in_bvhTree, const float4* in_verts, unsigned int in_verts_count,
                       const float in_innerC = 1.0, const float in_leafC = 1.0);
}



