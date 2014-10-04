#include "bem_matlab.h"
#include "bem_common.hpp"

/*
 * Please note, that although we are implementing a C-interface function
 * this file is written in C++. This module is used by the mx_bem_dlp
 * module.
 */

#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#define TRIANGLE(i) \
    { \
      coordinates[ 0*nC + static_cast<int>(elements[ 0*nE + i ]) - 1 ], \
      coordinates[ 1*nC + static_cast<int>(elements[ 0*nE + i ]) - 1 ], \
      coordinates[ 2*nC + static_cast<int>(elements[ 0*nE + i ]) - 1 ], \
      coordinates[ 0*nC + static_cast<int>(elements[ 1*nE + i ]) - 1 ], \
      coordinates[ 1*nC + static_cast<int>(elements[ 1*nE + i ]) - 1 ], \
      coordinates[ 2*nC + static_cast<int>(elements[ 1*nE + i ]) - 1 ], \
      coordinates[ 0*nC + static_cast<int>(elements[ 2*nE + i ]) - 1 ], \
      coordinates[ 1*nC + static_cast<int>(elements[ 2*nE + i ]) - 1 ], \
      coordinates[ 2*nC + static_cast<int>(elements[ 2*nE + i ]) - 1 ], \
    };

#ifdef __cplusplus
extern "C"
{
#endif

void
build_dlp_for_C(
  int nC,
  const double* coordinates,
  int nE,
  const double* elements,
  const double* areas,
  double minus_identity_factor,
  double* dlp)
{
  const std::vector<double> weights =
    ::millg::bem_util::get_gauss_weights();

  for (int i = 0; i < nE*nC; ++i)
  {
    dlp[i] = 0.;
  }

  for (int i = 0; i < nE; ++i)
  {
    const double triangle[9] = TRIANGLE(i);

    std::vector<Eigen::Vector3d> evaluation_points =
      ::millg::bem_util::compute_gauss_points(
          triangle, triangle + 3, triangle + 6);

    for (int j = 0; j < nE; ++j)
    {
      const double inner[9] = TRIANGLE(j);

      for (int k = 0; k < 3; ++k)
      {
        int vertex_id = elements[ j + k*nE ] - 1;

        dlp[ vertex_id * nE + i ] +=
          2 * areas[i] *
          ::millg::dlp_compute_local_contributions_for_patch_element(
            inner + 3*k,
            inner,
            evaluation_points,
            weights);

        if (i == j)
        {
          dlp[ vertex_id * nE + i ] +=
            minus_identity_factor * areas[j] / 3.;
        }
      }
    }
  }
}

#ifdef __cplusplus
}
#endif

