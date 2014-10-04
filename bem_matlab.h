#ifndef MILLG_BEM_MATLAB_HPP_GUARD_
#define MILLG_BEM_MATLAB_HPP_GUARD_

/**
 * @file bem_matlab.hpp
 * @brief A plain C interface to compute the double layer potential.
 */

/**
 * @brief Builds the double-layer potential for Laplace in three dimensions
 *        using a coordinates/elements mesh.
 *
 * @param nC The number of vertices in the mesh.
 * @param coordinates The coordinates of the vertices in the mesh, i.e. an
 *          3*nC array. The coordinates of the first vertex are located in
 *          coordinates[0], coordinates[nC] and coordinates[2*nC].
 * @param nE The number of elements (faces) in the mesh.
 * @param elements An 3*nE array that describes the elements in the mesh.
 *          The entries describing the first face are located in
 *          elements[0], elements[nE] and elements[2*nE].
 * @param areas An nE array containing the areas of the individual elements.
 * @param dlp An nE x nC array that contains the double layer potential matrix.
 *          The matrix is stored column-wise.
 *
 * @returns Fills the dlp array.
 */

#ifdef __cplusplus
extern "C"
{
#endif

extern void
build_dlp_for_C(
  int nC,
  const double* coordinates,
  int nE,
  const double* elements,
  const double* areas,
  double minus_identity_factor,
  double* dlp);

#ifdef __cplusplus
}
#endif

#endif

