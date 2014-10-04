#ifndef MILLG_BEM_COMMON_HPP_GUARD_
#define MILLG_BEM_COMMON_HPP_GUARD_

#include <vector>

#include <Eigen/Core>

/**
 * @brief This file contains common functions for computing the
 *        double-layer potential for the Laplace equation in
 *        three dimensions.
 *
 * This file is used by some other projects that do not use the same set
 * of dependencies as millg does. It contains all difficult parts to compute
 * the double-layer-potential for Laplace equation in three dimensions, but
 * leaves out some very mesh-specific stuff as well as the actual function
 * that assembles the matrix.
 *
 * For examples on how to use this module, you may want to take a look at
 * either the bem module or the bem_matlab module.
 */

namespace millg
{
  /**
   * @brief Contains functions useful for implementing the double-layer
   *        potential.
   */
  namespace bem_util
  {
    /**
     * @brief Computes transformed Gauss points for a given triangle.
     *
     * @param a The coordinates of the first point of the triangle.
     * @param b The coordinates of the second point of the triangle.
     * @param c The coordinates of the third point of the triangle.
     *
     * @returns We use a 7-point Gauss quadrature rule. This function
     *          transforms some quadrature points that are given on a
     *          reference triangle to quadrature points on the triangle
     *          given by (a,b,c).
     */
    std::vector<Eigen::Vector3d>
    compute_gauss_points(
      const double a[3], const double b[3], const double c[3]);

    /**
     * @brief Returns the Gauss weights for the same quadrature rule that
     *        compute_gauss_points returns the evaluation points for.
     */
    std::vector<double>
    get_gauss_weights();

    /**
     * @brief Helper method for quadrature.
     *
     * @param function An object \f$f\f$ that provides the operation
     *          TWeight operator()(const TEval&)
     * @param evaluation_points A list of evaluation points
     *          \f$(e_i)_{i=1}^n\f$.
     * @param weights A list of summation weights \f$(w_i)_{i=1}^n\f$
     *
     * The vectors evaluation_points and weights must have the same size.
     *
     * @returns The sum \f$\sum_{i=1}^n w_i \cdot f(e_i)\f$.
     */
    template < typename TEval, typename TWeight, typename TFunc >
    double quad(
      const TFunc& function,
      const std::vector<TEval>& evaluation_points,
      const std::vector<TWeight>& weights);
  }

  /**
   * @param boundary_vertex The coordinates of a boundary vertex \f$j\f$.
   * @param coordinates A given triangle \f$\tau\f$ that is oriented in the
   *        same way as all other triangles in the boundary mesh.
   * @param evaluation_points A list of evaluation points
   *          \f$(e_i)_{i=1}^n\f$.
   * @param weights A list of summation weights \f$(w_i)_{i=1}^n\f$
   *
   * The vectors evaluation_points and weights must have the same size.
   *
   * @returns The sum
   *  \f[
          \sum_{i=1}^n w_i \int_\tau \frac{x_i - y, n(y))}{|x_i-y|^3}
            \phi_j(y) dy
      \f]
   *  where \f$n(y)\f$ denotes the outer normal vector at \f$y\f$ and
   *  \f$\phi_j\f$ denotes the hat-function at the \f$j\f$-th node.
   */
  double
  dlp_compute_local_contributions_for_patch_element(
    const double boundary_vertex[3],
    const double coordinates[9],
    const std::vector<Eigen::Vector3d>& evaluation_points,
    const std::vector<double>& weights);
}

#endif

