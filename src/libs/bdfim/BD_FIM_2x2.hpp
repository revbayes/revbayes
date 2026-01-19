#ifndef BD_FIM_2X2
#define BD_FIM_2X2

#include <array>
#include <cmath>
#include <iostream>

// Detect Boost Gauss–Kronrod availability
#if defined(__has_include)
#  if __has_include(<boost/math/quadrature/gauss_kronrod.hpp>)
#    include <boost/math/quadrature/gauss_kronrod.hpp>
#    define BDGRADS_HAS_BOOST 1
#  else
#    define BDGRADS_HAS_BOOST 0
#  endif
#else
#  define BDGRADS_HAS_BOOST 0
#endif

namespace BDGrads {

  namespace TwoByTwo {

    //////////////////////////////////
    // helper structs and functions //
    //////////////////////////////////

    struct PComponents {
      double P;
      double dP[2];     // [dP/dλ, dP/dμ]
      double d2P[2][2]; // second derivatives
    };

    inline void computeP(const double& lambda, const double& mu, const double& rho, const double& t, const double& A, const double& C, PComponents& out) {

      double Et = std::exp(-(lambda - mu) * t);
      
      // B and derivatives (same as in your R function)
      double B   = rho * lambda + C * Et;
      
      double Bl  = rho + (1.0 - rho) * Et - t * C * Et;
      double Bm  = -Et + t * C * Et;
      
      double Bll = -2.0 * t * (1.0 - rho) * Et + t * t * C * Et;
      double Bmm = t * Et * (t * C - 2.0);
      double Blm = t * Et * (2.0 - rho - t * C);
      
      // P = A/B
      double P = A / B;
      
      // First derivatives: P_i = (A_i B - A B_i)/B^2
      double Al = rho;
      double Am = -rho;
      
      double Pl = (Al * B - A * Bl) / (B * B);
      double Pm = (Am * B - A * Bm) / (B * B);
      
      // Helper for second derivatives
      // P_ij = (Aij*B + Ai*Bj - Aj*Bi - A*Bij)/B^2 - 2*((Ai*B - A*Bi)*Bj)/B^3
      double A_ll = 0.0, A_mm = 0.0, A_lm = 0.0;
      
      double Pij_ll, Pij_mm, Pij_lm;
      
      // (λ,λ)
      {
          double num = A_ll*B + Al*Bl - Al*Bl - A*Bll;
          double term2 = (Al*B - A*Bl)*Bl;
          Pij_ll = num/(B*B) - 2.0*term2/(B*B*B);
      }
      
      // (μ,μ)
      {
          double num = A_mm*B + Am*Bm - Am*Bm - A*Bmm;
          double term2 = (Am*B - A*Bm)*Bm;
          Pij_mm = num/(B*B) - 2.0*term2/(B*B*B);
      }
      
      // (λ,μ)
      {
          double num = A_lm*B + Al*Bm - Am*Bl - A*Blm;
          double term2 = (Al*B - A*Bl)*Bm;
          Pij_lm = num/(B*B) - 2.0*term2/(B*B*B);
      }
      
      out.P      = P;
      out.dP[0]  = Pl;
      out.dP[1]  = Pm;
      
      out.d2P[0][0] = Pij_ll;
      out.d2P[1][1] = Pij_mm;
      out.d2P[0][1] = Pij_lm;
      out.d2P[1][0] = Pij_lm;

    }

    struct FIMIntegrand {

        double lambda;
        double mu;
        double rho;
        double A;
        double C;
        double vT;

        FIMIntegrand(const double& lambda_, const double& mu_, const double& rho_, const double& A_, const double& C_, const double& vT_) : lambda(lambda_), mu(mu_), rho(rho_), A(A_), C(C_), vT(vT_) {}
        
        void operator()(double t, double out[3]) const {
            
            // compute derivatives of P
            PComponents Pc;
            computeP(lambda, mu, rho, t, A, C, Pc);

            // retrieve calculated values
            const double  &P          = Pc.P;
            const double (&d2P)[2][2] = Pc.d2P;
            const double  &Pl         = Pc.dP[0];
            const double  &Pm         = Pc.dP[1];
            
            // p1(t) = (1/rho) * P^2 * exp(-(lambda-mu)t)
            double Et = std::exp(-(lambda - mu) * t);
            double p1 = (1.0 / rho) * P * P * Et;
            
            // g(t) = lambda * p1 / vT
            double g = lambda * p1 / vT;
            
            // Hessian of log p1(t):
            // H = 2 * (P2/P - dP dP^T / P^2)
            double H[2][2];
            H[0][0] = 2.0 * (d2P[0][0] / P - (Pl * Pl) / (P * P));
            H[1][1] = 2.0 * (d2P[1][1] / P - (Pm * Pm) / (P * P));
            H[0][1] = 2.0 * (d2P[0][1] / P - (Pl * Pm) / (P * P));
            H[1][0] = H[0][1];
            
            out[0] = H[0][0] * g; // λλ
            out[1] = H[0][1] * g; // λμ = μλ
            out[2] = H[1][1] * g; // μμ

        }

    };

    //////////////////////////
    // integraton functions //
    //////////////////////////

    inline void eval_f(const FIMIntegrand& f, double x, double out[3]) {
        f(x, out);
    }

    inline void adaptive_simpson_rec(const FIMIntegrand& f,
                              double a, double b,
                              const double fa[3],
                              const double fm[3],
                              const double fb[3],
                              const double S[3],
                              double tol,
                              int depth,
                              int maxDepth,
                              double result[3]) {

        double m  = 0.5 * (a + b);
        double h  = 0.5 * (b - a);
        double lm = 0.5 * (a + m);
        double rm = 0.5 * (m + b);
        
        double flm[3], frm[3];
        eval_f(f, lm, flm);
        eval_f(f, rm, frm);
        
        double Sleft[3], Sright[3], S2[3], err[3];
        double maxErr = 0.0;
        
        for (int k = 0; k < 3; ++k) {
            Sleft[k]  = (h / 6.0) * (fa[k] + 4.0 * flm[k] + fm[k]);
            Sright[k] = (h / 6.0) * (fm[k] + 4.0 * frm[k] + fb[k]);
            S2[k]     = Sleft[k] + Sright[k];
            err[k]    = std::fabs(S2[k] - S[k]);
            if (err[k] > maxErr) maxErr = err[k];
        }
        
        if (depth >= maxDepth || maxErr < 15.0 * tol) {
            // Accept with Richardson extrapolation
            for (int k = 0; k < 3; ++k) {
                result[k] = S2[k] + (S2[k] - S[k]) / 15.0;
            }
            return;
        } else {
            double S_left_res[3], S_right_res[3];
            
            adaptive_simpson_rec(f, a, m,
                                 fa, flm, fm, Sleft,
                                 0.5 * tol, depth + 1, maxDepth,
                                 S_left_res);
            
            adaptive_simpson_rec(f, m, b,
                                 fm, frm, fb, Sright,
                                 0.5 * tol, depth + 1, maxDepth,
                                 S_right_res);
            
            for (int k = 0; k < 3; ++k) {
                result[k] = S_left_res[k] + S_right_res[k];
            }
            return;
        }
    }

    // Wrapper: integrate f over [a,b] with tolerance tol
    inline void adaptive_simpson(const FIMIntegrand& f,
                          double a, double b,
                          double tol,
                          int maxDepth,
                          double result[3]) {
        double fa[3], fm[3], fb[3], S[3];
        
        eval_f(f, a, fa);
        eval_f(f, b, fb);
        double m = 0.5 * (a + b);
        eval_f(f, m, fm);
        
        double h = b - a;
        for (int k = 0; k < 3; ++k) {
            S[k] = (h / 6.0) * (fa[k] + 4.0 * fm[k] + fb[k]);
        }
        
        adaptive_simpson_rec(f, a, b,
                             fa, fm, fb, S,
                             tol, 0, maxDepth,
                             result);

    }

    /////////////////////////////////
    // integration backend switch  //
    /////////////////////////////////

    inline void integrate_1d(const FIMIntegrand& f,
                             double a, double b,
                             double tol,
                             int maxDepth,
                             double result[3],
                             bool force_simpson) {

#if BDGRADS_HAS_BOOST

      if (force_simpson) {

        // fallback: use original adaptive Simpson integrator
        adaptive_simpson(f, a, b, tol, maxDepth, result);

      } else {

        using boost::math::quadrature::gauss_kronrod;

        // integrate each component separately as a scalar function
        auto f11 = [&](double t){ double out[3]; f(t, out); return out[0]; };
        auto f12 = [&](double t){ double out[3]; f(t, out); return out[1]; };
        auto f22 = [&](double t){ double out[3]; f(t, out); return out[2]; };

        result[0] = gauss_kronrod<double, 15>::integrate(f11, a, b, tol);
        result[1] = gauss_kronrod<double, 15>::integrate(f12, a, b, tol);
        result[2] = gauss_kronrod<double, 15>::integrate(f22, a, b, tol);

        (void)maxDepth; // unused in this branch

      }


#else

      // fallback: use original adaptive Simpson integrator
      adaptive_simpson(f, a, b, tol, maxDepth, result);

#endif

    }

    ///////////////////////////////////
    // things that actually get used //
    ///////////////////////////////////

    inline std::array<std::array<double, 2>, 2> BD_FIM(double lambda, double mu, double rho, double T, double tol = 1e-6, int maxDepth = 15, bool force_simpson = false) {

      // precompute
      double A = rho * (lambda - mu);
      double C = lambda * (1.0 - rho) - mu;

      // get derivatives of Q etc.
      PComponents QT;
      computeP(lambda, mu, rho, T, A, C, QT);
      const double  &Q          = QT.P;
      const double (&dQ)[2]     = QT.dP;
      const double (&d2Q)[2][2] = QT.d2P;
      
      double FT = std::exp(-(lambda - mu) * T);
      double vT = 1.0 - Q * FT / rho;

      // first derivatives of v
      double v1[2];
      v1[0] = -(FT / rho) * (dQ[0] - T * Q); // dvT/dλ
      v1[1] = -(FT / rho) * (dQ[1] + T * Q); // dvT/dμ

      // second derivatives of v
      double v2[2][2];
      v2[0][0] = -(FT / rho) * (d2Q[0][0] - 2.0 * T * dQ[0] + T * T * Q);
      v2[1][1] = -(FT / rho) * (d2Q[1][1] + 2.0 * T * dQ[1] + T * T * Q);
      v2[0][1] = v2[1][0] = -(FT / rho) * (d2Q[0][1] + T * dQ[0] - T * dQ[1] - T * T * Q);

      // fixed terms
      double denom = 1.0 - vT;
      double denom2 = denom * denom;
      double H_log1mv[2][2];
      H_log1mv[0][0] = -(v2[0][0]) / denom - (v1[0] * v1[0]) / denom2;
      H_log1mv[0][1] = -(v2[0][1]) / denom - (v1[0] * v1[1]) / denom2;
      H_log1mv[1][0] = -(v2[1][0]) / denom - (v1[1] * v1[0]) / denom2;
      H_log1mv[1][1] = -(v2[1][1]) / denom - (v1[1] * v1[1]) / denom2;

      double H_loglambda[2][2];
      H_loglambda[0][0] = -1.0 / (lambda * lambda);
      H_loglambda[0][1] = 0.0;
      H_loglambda[1][0] = 0.0;
      H_loglambda[1][1] = 0.0;

      double En1 = vT / (1.0 - vT); // the expected number of tips

      // compute integrals (via Boost if available, otherwise Simpson)
      FIMIntegrand f(lambda, mu, rho, A, C, vT);
      double integral[3];
      integrate_1d(f, 0.0, T, tol, maxDepth, integral, force_simpson);

      // compute the information matrix
      std::array<std::array<double, 2>, 2> fim{{}};
      fim[0][0] = -(H_log1mv[0][0] + En1 * (H_loglambda[0][0] + integral[0]));
      fim[0][1] = fim[1][0] = -(H_log1mv[0][1] + En1 * (H_loglambda[0][1] + integral[1]));
      fim[1][1] = -(H_log1mv[1][1] + En1 * (H_loglambda[1][1] + integral[2]));

      return fim;

    }

  }

};

#endif
