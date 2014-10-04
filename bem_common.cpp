#include "bem_common.hpp"

#include <Eigen/Dense>

namespace millg
{
  namespace private__
  {
    const double bem_eps = 1e-13;
  
    double dlp_v(
      double alpha,
      double s, double p, double q)
    {
      assert(alpha == alpha);
      assert(s == s);
      assert(p == p);
      assert(q == q);
      if (std::fabs(s-p) < bem_eps)
        return 0.;
        
      double result = (std::sqrt( (1+alpha*alpha)*(s-p)*(s-p)+q*q ) - q) /
                    (std::sqrt(1+alpha*alpha) * (s-p));
      assert(result == result);
      return result;
    }
  
    /**
     * Computes the integral
     *
     * \f[
          \frac{2 u_x}{u_x^2+\alpha^2 q^2} \int \frac{D v + E}{v^2 + A v + B} dv,
       \f]
     * where
     * \f[
          A := -\frac{2\alpha \sqrt{1+\alpha^2}a}{u_x^2+\alpha^2 q^2}
                \left(\frac{t_x-\alpha s_x}{1+\alpha^2} + \sigma q\right), \\
          B := \frac{1+\alpha^2}{u_x^2+\alpha^2 q^2}
                \left( \frac{t_x-\alpha s_x}{1+\alpha^2} + \sigma q\right)^2, \\
          D := \sigma \frac12 \left( u_x^2+\alpha^2 a^2\right), \\
          E := \frac{1}{2\sqrt{1+\alpha^2}}
                  \left(s_\tau-s_x+\sigma \alpha q\right)
                  \left(t_x-\alpha s_x + \sigma\left(1+\alpha^2\right)q\right)
       \f]
     *
     * @param sgn Must be either 1 or -1.
     *
     * @returns The value as stated above.
     */
    double
    dlp_remaining_integral_helper(
      double alpha,
      double q,
      double u_x, double s_x, double t_x,
      double s_tau, int sgn,
      double v)
    {
      assert(sgn == 1 || sgn == -1);
  
      // Precondition: None of the input parameters is NaN.
      assert(alpha == alpha);
      assert(q == q);
      assert(u_x == u_x);
      assert(s_x == s_x);
      assert(t_x == t_x);
      assert(s_tau == s_tau);
      assert(v == v);
  
      double aux_1 = 1. / (u_x*u_x + alpha*alpha*q*q);
      double aux_2 = ((t_x - alpha*s_x) / (1 + alpha*alpha) + sgn*q);
  
      // Invariant: aux_1 and aux_2 must not be NaN.
      assert(aux_1 == aux_1);
      assert(aux_2 == aux_2);
  
      double A = -2 * alpha * q * std::sqrt(1+alpha*alpha) * aux_1 * aux_2;
      double B = (1+alpha*alpha) * aux_1 * aux_2 * aux_2;
      assert(4*B > A*A);
  
      double G = std::sqrt(B-.25*A*A);
      assert(G >= 0);
    
      if ( std::fabs( u_x ) > bem_eps )
      {
        double aux_atan = 0.;
        if ( std::fabs(G) < bem_eps )
        {
          if ( 2*v+A > bem_eps )
            aux_atan = +M_PI / 2.;
          else if ( 2*v+A < bem_eps )
            aux_atan = -M_PI / 2.;
          else
            assert(0);
        }
        else
        {
          aux_atan = std::atan( (2*v+A) / (2*G) );
        }
  
        // TODO: Replace by something more efficient:
        double u_x_sgn = u_x / std::fabs(u_x);
  
        assert(v*v+A*v+B > bem_eps);
        return -sgn * .5 * u_x * std::log(v*v+A*v+B) +
                sgn * (s_tau - s_x) * u_x_sgn * aux_atan;
      }
      else
      {
        return 0.;
      }
    }
  
    double
    dlp_remaining_integral_helper(
      double alpha,
      double s, double p, double q,
      double u_x, double s_x, double t_x,
      double s_tau, int sgn)
    {
      return dlp_remaining_integral_helper(alpha,q,u_x,s_x,t_x,s_tau,sgn,
        dlp_v(alpha,s,p,q));
    }
  
    /**
     * Computes the integral
     *
     * \f[
         I := \int \frac{
            \left(\alpha \left(s_\tau-s_x\right) - \left(\alpha s_x-t_x\right)\right)
              \left(s-s_x\right) +
            \left(s_tau-s_x\right)\left(\alpha s_x-t_x\right) + 
            \alpha u_x^2}
          {\left(\left(s-s_x\right)^2+u_x^2\right)
            \sqrt{\left(1+\alpha^2\right)\left(s-p\right)^2+q^2}} ds
       \f]
     */
    double
    dlp_remaining_integral(
      double alpha,
      double q,
      double u_x, double s_x, double t_x,
      double s_tau,
      double v)
    {
      double aux1 =
        dlp_remaining_integral_helper(alpha, q, u_x, s_x, t_x, s_tau,  1, v);
      double aux2 =
        dlp_remaining_integral_helper(alpha, q, u_x, s_x, t_x, s_tau, -1, v);
      assert(aux1 == aux1);
      assert(aux2 == aux2);
  
      double result = aux1 + aux2;
      assert(result == result);
  
      return result;
    }
  
    double
    dlp_remaining_integral(
      double alpha,
      double s, double p, double q,
      double u_x, double s_x, double t_x,
      double s_tau)
    {
      return dlp_remaining_integral(alpha, q, u_x, s_x, t_x, s_tau,
                                      dlp_v(alpha, s, p, q));
    }
  
    /* Computes
     *
     * \f[
          -\alpha u_x \int \frac{1}{\sqrt{(1+\alpha^2)(s-p)^2+q^2}} ds
       \f]
     *
     * using its anti-derivative.
     */
    double
    dlp_compute_easy_integral(
      double alpha,
      double s,
      double p,
      double q,
      double u_x)
    {
      assert( std::fabs(alpha) > 0 );
      assert( std::fabs(q) > 0 || (s-p) > 0 );
  
      double val = (u_x * alpha / std::sqrt(1+alpha*alpha) *
                   std::log( std::sqrt(1+alpha*alpha) * (s-p)
                   + std::sqrt((1+alpha*alpha)*(s-p)*(s-p)+q*q) ));
  
      assert( val == val );
      return val;
    }
  
    /* Computes
     *
     * \f[
          u_x I -
          u_x \frac{\alpha}{\sqrt{1+\alpha^2}} \log\left(
            \sqrt{1+\alpha^2}(s-p) + \sqrt{(1+\alpha^2)(s-p)^2+q^2}\right)
       \f]
     */
    double
    dlp_compute_F(
      double alpha,
      double s,
      double u_x, double s_x, double t_x,
      double s_tau)
    {
      if ( std::fabs( u_x ) < bem_eps )
      {
        return 0.;
      }
      else
      {
        assert(alpha == alpha);
        assert(s == s);
        assert(u_x == u_x);
        assert(s_x == s_x);
        assert(t_x == t_x);
        assert(s_tau == s_tau);
  
        double p = (alpha * t_x + s_x) / ( 1+alpha*alpha );
        double q = std::sqrt(u_x*u_x +
                    (t_x-alpha*s_x)*(t_x-alpha*s_x) / (1+alpha*alpha));
  
        assert(p == p);
        assert(q == q);
        assert( !(std::fabs(q) < bem_eps) || (std::fabs(u_x) < bem_eps) );
  
        double aux1 = 0.;
        if ( std::fabs( alpha ) > bem_eps && std::fabs( q ) > bem_eps )
        {
          aux1 = dlp_compute_easy_integral(alpha, s, p, q, u_x);
        }
    
        double remaining_integral =
          dlp_remaining_integral(alpha, s, p, q, u_x, s_x, t_x, s_tau);
        assert(remaining_integral == remaining_integral);
  
        return remaining_integral - aux1;
      }
    }
  
    /**
     * Computes the integral
     *
     * \f[
        \frac{1}{4\pi} \int_\tau \frac{\left(x-y, n(y)\right)}{|x-y|^3}
          \psi_{\tau,i}(y) ds_y
       \f]
     */
    double
    dlp_compute_local_contribution(
      const double triangle_coord[9],
      const Eigen::Vector3d& x)
    {
      Eigen::Vector3d A(triangle_coord[0],triangle_coord[1],triangle_coord[2]);
      Eigen::Vector3d B(triangle_coord[3],triangle_coord[4],triangle_coord[5]);
      Eigen::Vector3d C(triangle_coord[6],triangle_coord[7],triangle_coord[8]);
  
      Eigen::Vector3d r2 = C - B;
  
      double t_tau = r2.norm();
      assert(t_tau == t_tau);
      assert(t_tau > 0.);
  
      r2 = r2 / t_tau;
  
      double t_star = ( A - B ).dot(r2);
      assert(t_star == t_star);
  
      Eigen::Vector3d x_star = B + t_star * r2;
  
      Eigen::Vector3d r1 = x_star - A;
  
      double s_tau = r1.norm();
      assert(s_tau == s_tau);
      assert(s_tau > 0.);
  
      r1 = r1 / s_tau;
  
      Eigen::Vector3d n = r1.cross(r2);
  
      double alpha1 = -t_star / s_tau;
      double alpha2 = (t_tau - t_star) / s_tau;
      double s_x    = (x-A).dot(r1);
      double t_x    = (x-A).dot(r2);
      double u_x    = (x-A).dot(n);
  
      return .25 * (
                dlp_compute_F( alpha2, s_tau, u_x, s_x, t_x, s_tau )
              - dlp_compute_F( alpha2, 0.,    u_x, s_x, t_x, s_tau )
              - dlp_compute_F( alpha1, s_tau, u_x, s_x, t_x, s_tau )
              + dlp_compute_F( alpha1, 0.,    u_x, s_x, t_x, s_tau ))
                  / (M_PI * s_tau);
    }
  
    double
    dlp_compute_local_contribution_(
      double triangle_coord[9],
      double vertex[3])
    {
      Eigen::Vector3d v;
      v << vertex[0], vertex[1], vertex[2];
      return dlp_compute_local_contribution(triangle_coord, v);
    }

    /**
     * @brief When initialized with a triangle \f$\tau\f$, an instance
     *        corresponds to the function \f$D_i(\tau, \cdot)\f$.
     */
    class dlp_compute_local_contribution_for_fixed_triangle {
      public:
        dlp_compute_local_contribution_for_fixed_triangle(
            const double triangle[9])
          : triangle_()
        {
          for (int i = 0; i < 9; ++i)
          {
            triangle_[i] = triangle[i];
          }
        }
        double operator()(const Eigen::Vector3d& evaluation_point) const
        {
          return dlp_compute_local_contribution(triangle_, evaluation_point);
        }
    
      private:
        double triangle_[9];
    };

  }

  namespace bem_util
  {
    std::vector<double>
    get_gauss_weights()
    {
      std::vector<double> weights(7);
      weights[0] = 0.062969590272414;
      weights[1] = 0.062969590272414;
      weights[2] = 0.062969590272414;
      weights[3] = 0.066197076394253;
      weights[4] = 0.066197076394253;
      weights[5] = 0.066197076394253;
      weights[6] = 0.112500000000000; 
  
      return weights;
    }
  
    std::vector<Eigen::Vector3d>
    compute_gauss_points(
      const double a[3], const double b[3], const double c[3])
    {
      const double quadx[7] = {
          0.101286507323456,
          0.797426985353087,
          0.101286507323456,
          0.470142064105115,
          0.470142064105115,
          0.059715871789770,
          0.333333333333333,
      };
      const double quady[7] = {
          0.101286507323456,
          0.101286507323456,
          0.797426985353087,
          0.059715871789770,
          0.470142064105115,
          0.470142064105115,
          0.333333333333333,
      };
  
      Eigen::Vector3d A, B, C;
      A << a[0], a[1], a[2];
      B << b[0], b[1], b[2];
      C << c[0], c[1], c[2];
  
      Eigen::Vector3d x, y;
      x = B - A;
      y = C - A;
  
      std::vector<Eigen::Vector3d> quad_points(7);
      for (int i = 0; i < 7; ++i)
      {
        quad_points[i] = A + quadx[i]*x + quady[i]*y;
      }
  
      return quad_points;
    }

    template < typename TEval, typename TWeight, typename TFunc >
    double quad(
      const TFunc& function,
      const std::vector<TEval>& evaluation_points,
      const std::vector<TWeight>& weights)
    {
      size_t numQuad = evaluation_points.size();
      assert( numQuad == weights.size() );
    
      TWeight sum = 0.;
    
      for (size_t i = 0; i < numQuad; ++i)
      {
        sum += weights[i] * function(evaluation_points[i]);
      }
    
      return sum;
    }    
  } // End of namespace bem_util

  double
  dlp_compute_local_contributions_for_patch_element(
    const double boundary_vertex[3],
    const double coordinates[9],
    const std::vector<Eigen::Vector3d>& evaluation_points,
    const std::vector<double>& weights)
  {
    int runs = 0;

    double mycoordinates[9];
    for (int i = 0; i < 9; ++i)
      mycoordinates[i] = coordinates[i];
    
    while ((mycoordinates[0] != boundary_vertex[0]
            || mycoordinates[1] != boundary_vertex[1]
            || mycoordinates[2] != boundary_vertex[2])
            && runs < 5)
    {
      double tmp[3];
      for (int i = 0; i < 3; ++i)
      {
        tmp[i] = mycoordinates[i];
      }
      for (int i = 0; i < 6; ++i)
      {
        mycoordinates[i] = mycoordinates[3+i];
      }
      for (int i = 0; i < 3; ++i)
      {
        mycoordinates[6+i] = tmp[i];
      }
      runs += 1;
    }

    assert(runs < 5); // Consistency check: vertex is within the current
                      // patch element.

    // Quadrature:
    //   XXX: This minus does not belong here, but somewhere else!!!
    return -bem_util::quad(
              private__::dlp_compute_local_contribution_for_fixed_triangle(
                mycoordinates),
              evaluation_points, weights);
  }

} // End of namespace millg

