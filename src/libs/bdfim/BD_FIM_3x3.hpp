#ifndef BD_FIM_3X3
#define BD_FIM_3X3

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

  namespace ThreeByThree {

    //////////////////////////////////
    // helper structs and functions //
    //////////////////////////////////

    struct PComponents {
      double P;
      double dP[3];      // [dP/dλ, dP/dμ, dP/dρ]
      double d2P[3][3];  // second derivatives
    };

    inline void computeP(const double& lambda,
                         const double& mu,
                         const double& rho,
                         const double& t,
                         const double& A,
                         const double& C,
                         PComponents& out) {

      double Et = std::exp(-(lambda - mu) * t);

      // B and derivatives
      double B  = rho * lambda + C * Et;

      double Bl = rho + (1.0 - rho) * Et - t * C * Et;
      double Bm = -Et + t * C * Et;
      double Br = lambda * (1.0 - Et);

      double Bll = -2.0 * t * (1.0 - rho) * Et + t * t * C * Et;
      double Bmm = t * Et * (t * C - 2.0);
      double Blm = t * Et * (2.0 - rho - t * C);
      double Blr = 1.0 + (lambda * t - 1.0) * Et;
      double Bmr = -lambda * t * Et;
      double Brr = 0.0;

      // A and derivatives
      double Al = rho;
      double Am = -rho;
      double Ar = lambda - mu;

      double All = 0.0, Amm = 0.0, Arr = 0.0;
      double Alm = 0.0;
      double Alr = 1.0;
      double Amr = -1.0;

      double P = A / B;

      double Ai[3] = {Al, Am, Ar};
      double Bi[3] = {Bl, Bm, Br};

      double Aij[3][3] = {
        {All, Alm, Alr},
        {Alm, Amm, Amr},
        {Alr, Amr, Arr}
      };

      double Bij[3][3] = {
        {Bll, Blm, Blr},
        {Blm, Bmm, Bmr},
        {Blr, Bmr, Brr}
      };

      double B2 = B * B;
      double B3 = B2 * B;

      // First derivatives
      for (int i = 0; i < 3; ++i)
        out.dP[i] = (Ai[i] * B - A * Bi[i]) / B2;

      // Second derivatives
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          double num  = Aij[i][j] * B + Ai[i] * Bi[j] - Ai[j] * Bi[i] - A * Bij[i][j];
          double term2 = (Ai[i] * B - A * Bi[i]) * Bi[j];
          out.d2P[i][j] = num / B2 - 2.0 * term2 / B3;
        }
      }

      out.P = P;
    }

    //////////////////////////////////////////
    // integrand of the Fisher information //
    //////////////////////////////////////////

    struct FIMIntegrand {

      double lambda, mu, rho, A, C, vT;

      FIMIntegrand(double lambda_, double mu_, double rho_,
                   double A_, double C_, double vT_)
       : lambda(lambda_), mu(mu_), rho(rho_), A(A_), C(C_), vT(vT_) {}

      // out = [λλ, λμ, λρ, μμ, μρ, ρρ] contribution at time t
      void operator()(double t, double out[6]) const {

        PComponents Pc;
        computeP(lambda, mu, rho, t, A, C, Pc);

        const double &P = Pc.P;
        const double (&dP)[3] = Pc.dP;
        const double (&d2P)[3][3] = Pc.d2P;

        double Et = std::exp(-(lambda - mu) * t);
        double p1 = (1.0 / rho) * P * P * Et;
        double g  = lambda * p1 / vT;

        double invP  = 1.0 / P;
        double invP2 = invP * invP;

        double H[3][3];
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            H[i][j] = 2.0 * (d2P[i][j] * invP - dP[i] * dP[j] * invP2);

        // -log rho term → +1/ρ² in (ρ,ρ)
        H[2][2] += 1.0 / (rho * rho);

        out[0] = H[0][0] * g; // λλ
        out[1] = H[0][1] * g; // λμ
        out[2] = H[0][2] * g; // λρ
        out[3] = H[1][1] * g; // μμ
        out[4] = H[1][2] * g; // μρ
        out[5] = H[2][2] * g; // ρρ
      }
    };

    //////////////////////////
    // integration helpers  //
    //////////////////////////

    inline void eval_f(const FIMIntegrand& f, double x, double out[6]) {
      f(x, out);
    }

    inline void adaptive_simpson_rec(const FIMIntegrand& f,
                                     double a, double b,
                                     const double fa[6],
                                     const double fm[6],
                                     const double fb[6],
                                     const double S[6],
                                     double tol,
                                     int depth,
                                     int maxDepth,
                                     double result[6]) {

      double m  = 0.5 * (a + b);
      double h  = 0.5 * (b - a);
      double lm = 0.5 * (a + m);
      double rm = 0.5 * (m + b);

      double flm[6], frm[6];
      eval_f(f, lm, flm);
      eval_f(f, rm, frm);

      double Sleft[6], Sright[6], S2[6], err[6];
      double maxErr = 0.0;

      for (int k = 0; k < 6; ++k) {
        Sleft[k]  = (h / 6.0) * (fa[k] + 4.0 * flm[k] + fm[k]);
        Sright[k] = (h / 6.0) * (fm[k] + 4.0 * frm[k] + fb[k]);
        S2[k]     = Sleft[k] + Sright[k];
        err[k]    = std::fabs(S2[k] - S[k]);
        if (err[k] > maxErr) maxErr = err[k];
      }

      if (depth >= maxDepth || maxErr < 15.0 * tol) {
        for (int k = 0; k < 6; ++k)
          result[k] = S2[k] + (S2[k] - S[k]) / 15.0;
        return;
      }

      double left_res[6], right_res[6];

      adaptive_simpson_rec(f, a, m,
                           fa, flm, fm, Sleft,
                           0.5 * tol, depth+1, maxDepth,
                           left_res);

      adaptive_simpson_rec(f, m, b,
                           fm, frm, fb, Sright,
                           0.5 * tol, depth+1, maxDepth,
                           right_res);

      for (int k = 0; k < 6; ++k)
        result[k] = left_res[k] + right_res[k];
    }

    inline void adaptive_simpson(const FIMIntegrand& f,
                                 double a, double b,
                                 double tol,
                                 int maxDepth,
                                 double result[6]) {

      double fa[6], fb[6], fm[6], S[6];
      eval_f(f, a, fa);
      eval_f(f, b, fb);
      double m = 0.5 * (a + b);
      eval_f(f, m, fm);

      double h = b - a;
      for (int k = 0; k < 6; ++k)
        S[k] = (h / 6.0) * (fa[k] + 4.0 * fm[k] + fb[k]);

      adaptive_simpson_rec(f, a, b, fa, fm, fb, S, tol, 0, maxDepth, result);
    }

    inline void integrate_1d(const FIMIntegrand& f,
                             double a, double b,
                             double tol,
                             int maxDepth,
                             double result[6],
                             bool force_simpson) {

#if BDGRADS_HAS_BOOST

      if (force_simpson) {

        adaptive_simpson(f, a, b, tol, maxDepth, result);

      } else {

        using boost::math::quadrature::gauss_kronrod;

        auto f11 = [&](double t){ double out[6]; f(t, out); return out[0]; };
        auto f12 = [&](double t){ double out[6]; f(t, out); return out[1]; };
        auto f13 = [&](double t){ double out[6]; f(t, out); return out[2]; };
        auto f22 = [&](double t){ double out[6]; f(t, out); return out[3]; };
        auto f23 = [&](double t){ double out[6]; f(t, out); return out[4]; };
        auto f33 = [&](double t){ double out[6]; f(t, out); return out[5]; };

        result[0] = gauss_kronrod<double, 15>::integrate(f11, a, b, tol);
        result[1] = gauss_kronrod<double, 15>::integrate(f12, a, b, tol);
        result[2] = gauss_kronrod<double, 15>::integrate(f13, a, b, tol);
        result[3] = gauss_kronrod<double, 15>::integrate(f22, a, b, tol);
        result[4] = gauss_kronrod<double, 15>::integrate(f23, a, b, tol);
        result[5] = gauss_kronrod<double, 15>::integrate(f33, a, b, tol);

        (void)maxDepth;
      }

#else

      adaptive_simpson(f, a, b, tol, maxDepth, result);

#endif
    }

    ///////////////////////////////////
    // full 3×3 Fisher computation   //
    ///////////////////////////////////

    std::array<std::array<double, 3>, 3>
    inline BD_FIM(double lambda, double mu, double rho, double T,
           double tol = 1e-6, int maxDepth = 15, bool force_simpson = false) {

      double A = rho * (lambda - mu);
      double C = lambda * (1.0 - rho) - mu;

      PComponents QT;
      computeP(lambda, mu, rho, T, A, C, QT);

      const double &Q = QT.P;
      const double (&dQ)[3] = QT.dP;
      const double (&d2Q)[3][3] = QT.d2P;

      double FT = std::exp(-(lambda - mu) * T);

      double Ft1[3] = {-T * FT, T * FT, 0.0};
      double Ft2[3][3] =
        {{ T*T*FT,  -T*T*FT, 0.0 },
         { -T*T*FT,  T*T*FT, 0.0 },
         { 0.0,       0.0,    0.0 }};

      double s = Q * FT;

      double s1[3], s2[3][3];
      for (int i = 0; i < 3; ++i)
        s1[i] = dQ[i] * FT + Q * Ft1[i];

      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          s2[i][j] =
            d2Q[i][j] * FT +
            dQ[i] * Ft1[j] +
            dQ[j] * Ft1[i] +
            Q * Ft2[i][j];

      double vT = 1.0 - s / rho;

      double v1[3];
      v1[0] = -s1[0] / rho;
      v1[1] = -s1[1] / rho;
      v1[2] = (s - rho * s1[2]) / (rho * rho);

      double v2[3][3];

      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
          v2[i][j] = -s2[i][j] / rho;

      v2[0][2] = v2[2][0] = (s1[0] - rho * s2[0][2]) / (rho * rho);
      v2[1][2] = v2[2][1] = (s1[1] - rho * s2[1][2]) / (rho * rho);

      v2[2][2] =
        (-rho*rho * s2[2][2] + 2.0 * rho * s1[2] - 2.0 * s)
        / (rho * rho * rho);

      double denom = 1.0 - vT;
      double denom2 = denom * denom;

      double H_log1mv[3][3];
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          H_log1mv[i][j] =
            -(v2[i][j]) / denom -
            (v1[i] * v1[j]) / denom2;

      double H_loglambda[3][3] = {{0.0}};
      H_loglambda[0][0] = -1.0 / (lambda * lambda);

      double En1 = vT / (1.0 - vT);

      FIMIntegrand f(lambda, mu, rho, A, C, vT);
      double Iint[6];
      integrate_1d(f, 0.0, T, tol, maxDepth, Iint, force_simpson);

      std::array<std::array<double, 3>, 3> fim{};

      fim[0][0] = -(H_log1mv[0][0] + En1*(H_loglambda[0][0] + Iint[0]));
      fim[0][1] = fim[1][0] = -(H_log1mv[0][1] + En1*(Iint[1]));
      fim[0][2] = fim[2][0] = -(H_log1mv[0][2] + En1*(Iint[2]));
      fim[1][1] = -(H_log1mv[1][1] + En1*(Iint[3]));
      fim[1][2] = fim[2][1] = -(H_log1mv[1][2] + En1*(Iint[4]));
      fim[2][2] = -(H_log1mv[2][2] + En1*(Iint[5]));

      return fim;
    }

    /////////////////
    // BD gradient //
    /////////////////
    
    // compute vT and its derivatives using your existing computeP()
    inline void compute_vT(double lambda, double mu, double rho,
                            double T, double A, double C,
                            double &vT, double dv[3]) {

      PComponents Pc;
      computeP(lambda, mu, rho, T, A, C, Pc);

      double Q = Pc.P;
      double Ql = Pc.dP[0];
      double Qm = Pc.dP[1];
      double Qr = Pc.dP[2];

      double F = std::exp(-(lambda - mu) * T);

      // vT = 1 - (QF)/rho
      vT = 1.0 - (Q * F) / rho;

      // ∂F/∂λ = -T F, ∂F/∂μ = +T F, ∂F/∂ρ = 0
      double Fl = -T * F;
      double Fm =  T * F;
      double Fr = 0.0;

      // s = QF
      double s1[3];
      s1[0] = Ql * F + Q * Fl;
      s1[1] = Qm * F + Q * Fm;
      s1[2] = Qr * F + Q * Fr;

      // dv/dθ
      dv[0] = -s1[0] / rho;                 // λ
      dv[1] = -s1[1] / rho;                 // μ
      dv[2] = (Q * F) / (rho * rho) - s1[2] / rho; // ρ

    }

    /////////////////////////////////////////////////
    // MAIN GRADIENT FUNCTION TO ADD
    /////////////////////////////////////////////////

    inline std::array<double,3> BD_Joint_Gradient(double lambda,
                                                  double mu,
                                                  double rho,
                                                  double T,
                                                  const std::vector<double> &tvec) {

      const int n = (int)tvec.size() + 1;

      // constants used repeatedly
      double A = rho * (lambda - mu);
      double C = lambda * (1.0 - rho) - mu;

      // compute vT and derivatives
      double vT, dv[3];
      compute_vT(lambda, mu, rho, T, A, C, vT, dv);

      // Initialize gradient components
      double g_lambda = 0.0;
      double g_mu     = 0.0;
      double g_rho    = 0.0;

      // (1) contribution from log(1 - vT)
      g_lambda += -dv[0] / (1.0 - vT);
      g_mu     += -dv[1] / (1.0 - vT);
      g_rho    += -dv[2] / (1.0 - vT);

      // (2) contribution from (n - 1) log λ
      g_lambda += (n - 1) * (1.0 / lambda);

      // (3) contribution from sum_{i} log p1(t_i)
      for (double t : tvec) {
        PComponents Pc;
        computeP(lambda, mu, rho, t, A, C, Pc);

        double P = Pc.P;
        double Pl = Pc.dP[0];
        double Pm = Pc.dP[1];
        double Pr = Pc.dP[2];

        // log p1(t) = -log ρ + 2logP - (λ - μ)t
        g_lambda += 2.0 * (Pl / P) - t;
        g_mu     += 2.0 * (Pm / P) + t;
        g_rho    += -1.0/rho + 2.0 * (Pr / P);
      }

      return {g_lambda, g_mu, g_rho};

    }


  } // namespace ThreeByThree

} // namespace BDGrads

#endif