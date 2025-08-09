//
// This source code resides at www.jaeckel.org/LetsBeRational.7z .
//
// ======================================================================================
// Copyright © 2013-2024 Peter Jäckel.
//
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
//
// WARRANTY DISCLAIMER
// The Software is provided "as is" without warranty of any kind, either express or implied,
// including without limitation any implied warranties of condition, uninterrupted use,
// merchantability, fitness for a particular purpose, or non-infringement.
// ======================================================================================
//

#include "lets_be_rational.h"

// To cross-compile on a command line, you could just use something like
//
//   i686-w64-mingw32-g++ -w -fpermissive -shared -DNDEBUG -Ofast erf_cody.cpp lets_be_rational.cpp main.cpp normaldistribution.cpp rationalcubic.cpp XLFunctions.cpp XLOper.cpp -o LetsBeRational.xll -static-libstdc++ -static-libgcc -s
//
// To compile into a shared library on non-Windows systems, you can use
//
//   g++ -fPIC -shared -DNDEBUG -Ofast erf_cody.cpp lets_be_rational.cpp main.cpp normaldistribution.cpp rationalcubic.cpp XLFunctions.cpp XLOper.cpp -o LetsBeRational.so -s
//

#if defined(_MSC_VER)
# define NOMINMAX // to suppress MSVC's definitions of min() and max()
#endif

#include "normaldistribution.h"
#include "rationalcubic.h"
#include <float.h>
#include <cmath>
#include <algorithm>
#if defined(_WIN32) || defined(_WIN64)
# include <windows.h>
#endif

// Useful for RELEASE-only issues (old-fashioned pedestrian debugging).
//#define PRINTVAR( v ) { printf( #v ": %.17g\n",v); fflush(stdout); }

#if defined( LOG_BINARY_NESTING ) || defined( PRINTVAR ) || ( defined( SAFEGUARD_VEGA_AGAINST_ZERO_VOLATILITY ) && defined( LOG_VEGA_OVERFLOW ) )
#include <stdio.h>
#endif

#include <assert.h>

#include <tuple>

#define TWO_PI                            6.283185307179586476925286766559005768394338798750
#define SQRT_PI_OVER_TWO                  1.253314137315500251207882642405522626503493370305  // sqrt(pi/2) to avoid misinterpretation.
#define SQRT_THREE                        1.732050807568877293527446341505872366942805253810
#define SQRT_ONE_OVER_THREE               0.577350269189625764509148780501957455647601751270
#define TWO_PI_OVER_SQRT_TWENTY_SEVEN     1.209199576156145233729385505094770488189377498728 // 2*pi/sqrt(27)
#define SQRT_THREE_OVER_THIRD_ROOT_TWO_PI 0.938643487427383566075051356115075878414688769574 // √3 / ∛(2π)
#define PI_OVER_SIX                       0.523598775598298873077107230546583814032861566563

namespace {
  static const double SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);
  static const double FOURTH_ROOT_DBL_EPSILON = sqrt(SQRT_DBL_EPSILON);
  static const double EIGHTH_ROOT_DBL_EPSILON = sqrt(FOURTH_ROOT_DBL_EPSILON);
  static const double SIXTEENTH_ROOT_DBL_EPSILON = sqrt(EIGHTH_ROOT_DBL_EPSILON);
  static const double SQRT_DBL_MIN = sqrt(DBL_MIN);
  static const double SQRT_DBL_MAX = sqrt(DBL_MAX);

  // Define this to a positive non-denormalised number if you want to suppress positive results for for (positive) denormalised inputs.
  // NOT RECOMMENDED.
  //#define POSITIVE_DENORMALISATION_CUTOFF DBL_MIN

  static const double VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC = -DBL_MAX;
  static const double VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM = DBL_MAX;

#if defined( POSITIVE_DENORMALISATION_CUTOFF )
  inline bool is_below_horizon(double x) { return fabs(x) < POSITIVE_DENORMALISATION_CUTOFF; } // This weeds out denormalised (a.k.a. 'subnormal') numbers.
#endif

#if defined ( USE_ORIGINAL_REGION_DEFINITIONS )
# define ETA -10 /* η */
#else
# define ETA -13 /* η */
#endif
#define TAU (2 * SIXTEENTH_ROOT_DBL_EPSILON) /* τ */

#if defined( POSITIVE_DENORMALISATION_CUTOFF )
#define VOLATILITY_IS_CONSIDERED_ZERO(θx,s) (s <= fabs(θx) * POSITIVE_DENORMALISATION_CUTOFF)
#else
#define VOLATILITY_IS_CONSIDERED_ZERO(θx,s) (s <= 0)
#endif

#if defined ( USE_ORIGINAL_REGION_DEFINITIONS )
  // Denote h := θx/s and t := s/2. We evaluate the condition |h|>|η|, i.e., h<η  &&  t < τ+|h|-|η|  avoiding any divisions by s , where η = ETA  and τ = TAU .
#define IS_REGION_I(θx,s)   (θx < s * ETA && 0.5 * s * s + θx < s * (TAU + ETA))
#define IS_REGION_II(θx,s)  (0.5 * s < TAU)
#else
#define IS_REGION_I(θx,s)   (θx < s * ETA && s * (0.5 * s - (TAU + 0.5 + ETA)) + θx < 0) // h < η,  t < (τ+half) + (|h|-|η|)  ⇔  s/2 < (τ+half) - θ·x/s + η  ⇔  s·(s/2-(τ+half+η))+x < 0     using    h = θ·x/s and t = s/2
#define IS_REGION_II(θx,s)  (s * (s - (2 * TAU)) - θx / ETA < 0) // t < τ + (half/|η|)·|h|  ⇔  s < 2·τ + (θ·x/s)/η  ⇔  s·(s-2·τ)-θ·x/η < 0                  using    h = θ·x/s and t = s/2
#endif
// When b is more than, say, about 85% of bₘₐₓ = exp(θ·x/2), then b is dominated by the first of the two terms in the Black formula, and we retain more accuracy by not attempting to combine the two terms in any way.
// We evaluate the condition h+t>0.85  avoiding any division by s.
#define IS_REGION_III(θx,s) (s * (0.5 * s - 0.85) + θx > 0) // t+h > 0.85  ⇔  s/2+θ·x/s > 0.85  ⇔  s·(s/2-0.85)+x > 0      using    h = θ·x/s and t = s/2

#if defined( ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT ) || defined( ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER ) || defined( E )

  // See https://www.kernel.org/doc/Documentation/atomic_ops.txt for further details on this simplistic implementation of an atomic flag that is *not* volatile.
  typedef struct {
#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
    long data;
#else
    int data;
#endif
  } atomic_t;

#endif

#ifdef ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT

  static atomic_t implied_volatility_maximum_iterations = { 2 }; // (DBL_DIG*20)/3 ≈ 100 . Only needed when the iteration effectively alternates Householder/Halley/Newton steps and binary nesting due to roundoff truncation.

# define IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS implied_volatility_maximum_iterations.data

#else

# define IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS 2

#endif

#if defined( ENABLE_SAFE_GUARDING_IN_HOUSEHOLDER_EXPANSIONS ) // Tests show that there is in fact no need for this
  inline double householder3_factor(double ν, double h₂, double h₃) {
    // Safeguarding against out-of-bounds behaviour comprised by a change in sign with fallback to either Halley or Newton, whichever is admissible.
    // The Halley method is ν / (1 + eta) with eta := 0.5 · h₂ · ν. It should have the same sign as ν = 'newton' by itself.
    // Hence, 1 + eta <= 0 is a clear indicator that the current guess is not inside the domain of attraction of the Halley method and we should fall back to the Newton method.
    // The Housholder(3) method is designed and intended as an improvement over the Halley method whence, if Halley is already failing the sign check, we do not even dare to look at the Housholder(3) method.
    const double eta = 0.5 * h₂ * ν;
    if (eta > -1) {
      // The Housholder(3) method is ν * (1 + eta) / (1 + zeta) with zeta := ν * (h₂ + h₃ * ν / 6), and it should also have the same sign as ν = 'newton' by itself.
      const double zeta = ν * (h₂ + h₃ * ν / 6);
      if (zeta > -1)
        return (1 + eta) / (1 + zeta);
      return 1 / (1 + eta);
    }
    return 1;
  }
#else
  inline double householder3_factor(double ν, double h₂, double h₃) { return (1 + 0.5 * h₂ * ν) / (1 + ν * (h₂ + h₃ * ν * (1./6.))); }
#endif

  inline double householder4_factor(double ν, double h₂, double h₃, double h₄) { return (1 + ν * (h₂ + ν * h₃ * (1./6.))) / (1 + ν * (1.5 * h₂ + ν * (h₂ * h₂ * 0.25 + h₃ * (1./3.) + ν * h₄ * (1./24.)))); }

#ifdef ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER

  static atomic_t implied_volatility_maximum_householder_method_order = { 5 };

#if defined( ENABLE_SAFE_GUARDING_IN_HOUSEHOLDER_EXPANSIONS ) // Tests show that there is in fact no need for this
  inline double halley_factor(double ν, double h₂) {
    // Safeguarding against out-of-bounds behaviour comprised by a change in sign with fallback to Newton.
    // The Halley method is ν / (1 + eta) with eta := 0.5 · h₂ · ν. It should have the same sign as ν = 'newton' by itself.
    // Hence, 1 + eta <= 0 is a clear indicator that the current guess is not inside the domain of attraction of the Halley method and we should fall back to the Newton method.
    const double eta = 0.5 * h₂ * ν;
    if (eta > -1)
      return 1 / (1 + eta);
    return 1;
  }
#else
  inline double halley_factor(double ν, double h₂) { return 1 / (1 + 0.5 * h₂ * ν); }
#endif

  inline double householder_factor(double ν, double h₂, double h₃) {
    return implied_volatility_maximum_householder_method_order.data > 3 ? householder3_factor(ν, h₂, h₃) : (implied_volatility_maximum_householder_method_order.data > 2 ? halley_factor(ν, h₂) : 1);
  }
  inline double householder_factor(double ν, double h₂, double h₃, double h₄) {
    return implied_volatility_maximum_householder_method_order.data > 4 ? householder4_factor(ν, h₂, h₃, h₄) : householder_factor(ν, h₂, h₃);
  }

#else

  inline double householder_factor(double ν, double h₂, double h₃) { return householder3_factor(ν, h₂, h₃); }
  inline double householder_factor(double ν, double h₂, double h₃, double h₄) { return householder4_factor(ν, h₂, h₃, h₄); }

#endif

#ifdef ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT

  static atomic_t implied_volatility_output_type = { 0 };

  inline double implied_volatility_output(int count, double volatility) { return implied_volatility_output_type.data != 0 ? count : volatility; }

#else

  inline double implied_volatility_output(int /* count */, double volatility) { return volatility; }

#endif  

}

#ifdef ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT
int set_implied_volatility_maximum_iterations(int n) {
  if (n >= 0) {
#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
    long i = (long)n;
    InterlockedExchange(&(implied_volatility_maximum_iterations.data), i);
#elif defined( __x86__ ) || defined( __x86_64__ )
    implied_volatility_maximum_iterations.data = n;
#else
# error Atomic operations not implemented for this platform.
#endif
  }
  return (int)implied_volatility_maximum_iterations.data;
}
#endif

#ifdef ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER
int set_implied_volatility_householder_method_order(int order) {
  if (order >= 0) {
#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
    long i = (long)order;
    InterlockedExchange(&(implied_volatility_maximum_householder_method_order.data), i);
#elif defined( __x86__ ) || defined( __x86_64__ )
    implied_volatility_maximum_householder_method_order.data = order;
#else
# error Atomic operations not implemented for this platform.
#endif
  }
  return (int)implied_volatility_maximum_householder_method_order.data;
}
#endif  

#ifdef ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT
int set_implied_volatility_output_type(int type) {
  if (type >= 0) {
#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
    long i = (long)type;
    InterlockedExchange(&(implied_volatility_output_type.data), i);
#elif defined( __x86__ ) || defined( __x86_64__ )
    implied_volatility_output_type.data = type;
#else
# error Atomic operations not implemented for this platform.
#endif
  }
  return (int)implied_volatility_output_type.data;
}
#endif  

inline double normalised_intrinsic(double θx) {
#if defined( AVOID_SINH )
  if (θx <= 0)
    return 0;
  const double x2 = θx * θx;
  if (x2 < 98 * FOURTH_ROOT_DBL_EPSILON) // The factor 98 is computed from last coefficient: √√92897280 = 98.1749
    return θx * (1 + x2 * ((1.0 / 24.0) + x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2))));
  // Caution: using bₘₐₓ = exp(0.5 * θx) and 1/bₘₐₓ to reduce the number of exponentials below incurs a larger roundoff error.
  return exp(0.5 * θx) - exp(-0.5 * θx);
#else
  return θx <= 0 ? 0 : 2 * std::sinh(0.5 * θx);
#endif
}

inline double square(const double x) { return x * x; }

// Asymptotic expansion of the 'scaled normalised Black' function (in analogy to the 'scaled' complementary error function)
// 
//             bx  :=  b / ( ∂(b(x,s)/∂s )
// with
//
//              b   =  Φ(h+t)·exp(θ·x/2) - Φ(h-t)·exp(-θ·x/2)
// and
//              h   =  θ·x/s   and   t  =  s/2   for h < 0, |h| > 10, t ≲ |h| -9.8.
// which makes
//              b   =  Φ(h+t)·exp(h·t) - Φ(h-t)·exp(-h·t)
//
//                     exp(-(h²+t²)/2)
//                  =  ---------------  ·  [ Y(h+t) - Y(h-t) ] .
//                         √(2π)
// 
// Since the normalised vega is given by
// 
//      ∂(b(x,s)/∂s = exp(-(h²+t²)/2)/√(2π)
// 
// we thus have
//             bx  :=  Y(h+t) - Y(h-t)
// with
//            Y(z) := Φ(z)/φ(z)
//
// for large negative (t-|h|) by the aid of Abramowitz & Stegun (26.2.12) where Φ(z) = φ(z)/|z|·[1-1/z^2+...].
// We define
//                     r
//         A(h,t) :=  --- · [ Y(h+t) - Y(h-t) ]
//                     t
//
// with r := (h+t)·(h-t) and give an expansion for A(h,t) in q:=(h/r)² expressed in terms of e:=(t/h)² .
//
double asymptotic_expansion_of_scaled_normalised_black(double h, double t) {
  // h < η  &&  t < (τ+half) + (|h|-|η|)
  assert(h < -fabs(ETA) && h < TAU + 0.5 - h + ETA);
  // Note that e := (t/h)² ∈ (0,1).
  const double e = square(t / h), r = (h + t) * (h - t), q = square(h / r);
#if defined( USE_ORIGINAL_ASYMPTOTIC_EXPANSION )
  // 17th order asymptotic expansion of A(h,t) in q, sufficient for Φ(z) [and thus Y(z)] to have relative accuracy of 1.64E-16 for z <= η  with  η:=ETA.
  const double ω = (2.0 + q * (-6.0E0 - 2.0 * e + 3.0 * q * (1.0E1 + e * (2.0E1 + 2.0 * e) + 5.0 * q * (-1.4E1 + e * (-7.0E1 + e * (-4.2E1 - 2.0 * e)) + 7.0 * q * (1.8E1 + e * (1.68E2 + e * (2.52E2 + e * (7.2E1 + 2.0 * e))) + 9.0 * q * (-2.2E1 + e * (-3.3E2 + e * (-9.24E2 + e * (-6.6E2 + e * (-1.1E2 - 2.0 * e)))) + 1.1E1 * q * (2.6E1 + e * (5.72E2 + e * (2.574E3 + e * (3.432E3 + e * (1.43E3 + e * (1.56E2 + 2.0 * e))))) + 1.3E1 * q * (-3.0E1 + e * (-9.1E2 + e * (-6.006E3 + e * (-1.287E4 + e * (-1.001E4 + e * (-2.73E3 + e * (-2.1E2 - 2.0 * e)))))) + 1.5E1 * q * (3.4E1 + e * (1.36E3 + e * (1.2376E4 + e * (3.8896E4 + e * (4.862E4 + e * (2.4752E4 + e * (4.76E3 + e * (2.72E2 + 2.0 * e))))))) + 1.7E1 * q * (-3.8E1 + e * (-1.938E3 + e * (-2.3256E4 + e * (-1.00776E5 + e * (-1.84756E5 + e * (-1.51164E5 + e * (-5.4264E4 + e * (-7.752E3 + e * (-3.42E2 - 2.0 * e)))))))) + 1.9E1 * q * (4.2E1 + e * (2.66E3 + e * (4.0698E4 + e * (2.3256E5 + e * (5.8786E5 + e * (7.05432E5 + e * (4.0698E5 + e * (1.08528E5 + e * (1.197E4 + e * (4.2E2 + 2.0 * e))))))))) + 2.1E1 * q * (-4.6E1 + e * (-3.542E3 + e * (-6.7298E4 + e * (-4.90314E5 + e * (-1.63438E6 + e * (-2.704156E6 + e * (-2.288132E6 + e * (-9.80628E5 + e * (-2.01894E5 + e * (-1.771E4 + e * (-5.06E2 - 2.0 * e)))))))))) + 2.3E1 * q * (5.0E1 + e * (4.6E3 + e * (1.0626E5 + e * (9.614E5 + e * (4.08595E6 + e * (8.9148E6 + e * (1.04006E7 + e * (6.53752E6 + e * (2.16315E6 + e * (3.542E5 + e * (2.53E4 + e * (6.0E2 + 2.0 * e))))))))))) + 2.5E1 * q * (-5.4E1 + e * (-5.85E3 + e * (-1.6146E5 + e * (-1.77606E6 + e * (-9.37365E6 + e * (-2.607579E7 + e * (-4.01166E7 + e * (-3.476772E7 + e * (-1.687257E7 + e * (-4.44015E6 + e * (-5.9202E5 + e * (-3.51E4 + e * (-7.02E2 - 2.0 * e)))))))))))) + 2.7E1 * q * (5.8E1 + e * (7.308E3 + e * (2.3751E5 + e * (3.12156E6 + e * (2.003001E7 + e * (6.919458E7 + e * (1.3572783E8 + e * (1.5511752E8 + e * (1.0379187E8 + e * (4.006002E7 + e * (8.58429E6 + e * (9.5004E5 + e * (4.7502E4 + e * (8.12E2 + 2.0 * e))))))))))))) + 2.9E1 * q * (-6.2E1 + e * (-8.99E3 + e * (-3.39822E5 + e * (-5.25915E6 + e * (-4.032015E7 + e * (-1.6934463E8 + e * (-4.1250615E8 + e * (-6.0108039E8 + e * (-5.3036505E8 + e * (-2.8224105E8 + e * (-8.870433E7 + e * (-1.577745E7 + e * (-1.472562E6 + e * (-6.293E4 + e * (-9.3E2 - 2.0 * e)))))))))))))) + 3.1E1 * q * (6.6E1 + e * (1.0912E4 + e * (4.74672E5 + e * (8.544096E6 + e * (7.71342E7 + e * (3.8707344E8 + e * (1.14633288E9 + e * (2.07431664E9 + e * (2.33360622E9 + e * (1.6376184E9 + e * (7.0963464E8 + e * (1.8512208E8 + e * (2.7768312E7 + e * (2.215136E6 + e * (8.184E4 + e * (1.056E3 + 2.0 * e))))))))))))))) + 3.3E1 * (-7.0E1 + e * (-1.309E4 + e * (-6.49264E5 + e * (-1.344904E7 + e * (-1.4121492E8 + e * (-8.344518E8 + e * (-2.9526756E9 + e * (-6.49588632E9 + e * (-9.0751353E9 + e * (-8.1198579E9 + e * (-4.6399188E9 + e * (-1.6689036E9 + e * (-3.67158792E8 + e * (-4.707164E7 + e * (-3.24632E6 + e * (-1.0472E5 + e * (-1.19E3 - 2.0 * e))))))))))))))))) * q)))))))))))))))));
#else
#define A0   2
#define A1   -6-2*e
#define A2   30+e*(60+6*e)
#define A3   -2.1E2+e*(-1.05E3+e*(-6.3E2-30*e))
#define A4   1.89E3+e*(1.764E4+e*(2.646E4+e*(7.56E3+2.1E2*e)))
#define A5   -2.079E4+e*(-3.1185E5+e*(-8.7318E5+e*(-6.237E5+e*(-1.0395E5-1.89E3*e))))
#define A6   2.7027E5+e*(5.94594E6+e*(2.675673E7+e*(3.567564E7+e*(1.486485E7+e*(1.62162E6+2.079E4*e)))))
#define A7   -4.05405E6+e*(-1.2297285E8+e*(-8.1162081E8+e*(-1.73918745E9+e*(-1.35270135E9+e*(-3.6891855E8+e*(-2.837835E7-2.7027E5*e))))))
#define A8   6.891885E7+e*(2.756754E9+e*(2.50864614E10+e*(7.88431644E10+e*(9.85539555E10+e*(5.01729228E10+e*(9.648639E9+e*(5.513508E8+4.05405E6*e)))))))
#define A9   -1.30945815E9+e*(-6.678236565E10+e*(-8.013883878E11+e*(-3.4726830138E12+e*(-6.3665855253E12+e*(-5.2090245207E12+e*(-1.8699062382E12+e*(-2.671294626E11+e*(-1.178512335E10-6.891885E7*e))))))))
#define A10   2.749862115E10+e*(1.7415793395E12+e*(2.664616389435E13+e*(1.52263793682E14+e*(3.848890340295E14+e*(4.618668408354E14+e*(2.664616389435E14+e*(7.10564370516E13+e*(7.83710702775E12+e*(2.749862115E11+1.30945815E9*e)))))))))
#define A11   -6.3246828645E11+e*(-4.870005805665E13+e*(-9.2530110307635E14+e*(-6.74147946527055E15+e*(-2.24715982175685E16+e*(-3.71802806872497E16+e*(-3.14602375045959E16+e*(-1.34829589305411E16+e*(-2.77590330922905E15+e*(-2.4350029028325E14+e*(-6.95715115095E12-2.749862115E10*e))))))))))
#define A12   1.581170716125E13+e*(1.454677058835E15+e*(3.36030400590885E16+e*(3.04027505296515E17+e*(1.29211689751018875E18+e*(2.81916414002223E18+e*(3.289024830025935E18+e*(2.067387036016302E18+e*(6.8406188691715875E17+e*(1.12010133530295E17+e*(8.0007238235925E15+e*(1.89740485935E14+6.3246828645E11*e)))))))))))
#define A13   -4.2691609335375E14+e*(-4.624924344665625E16+e*(-1.2764791191277125E18+e*(-1.40412703104048375E19+e*(-7.41067044160255312E19+e*(-2.06151377739125569E20+e*(-3.17155965752500875E20+e*(-2.74868503652167425E20+e*(-1.33392067948845956E20+e*(-3.51031757760120938E19+e*(-4.6804234368016125E18+e*(-2.774954606799375E17+e*(-5.54990921359875E15-1.581170716125E13*e))))))))))))
#define A14   1.238056670725875E16+e*(1.5599514051146025E18+e*(5.06984206662245812E19+e*(6.66322100184665925E20+e*(4.27556680951827302E21+e*(1.47701398874267613E22+e*(2.89721974714909549E22+e*(3.31110828245610914E22+e*(2.2155209831140142E22+e*(8.55113361903654604E21+e*(1.83238577550783129E21+e*(2.02793682664898325E20+e*(1.01396841332449162E19+e*(1.733279339016225E17+4.2691609335375E14*e)))))))))))))
#define A15   -3.8379756792502125E17+e*(-5.56506473491280812E19+e*(-2.10359446979704147E21+e*(-3.25556286992399275E22+e*(-2.49593153360839444E23+e*(-1.04829124411552567E24+e*(-2.55352995361474201E24+e*(-3.72085793241005264E24+e*(-3.28310994036181115E24+e*(-1.74715207352587611E24+e*(-5.49104937393846778E23+e*(-9.76668860977197826E22+e*(-9.11557603578717971E21+e*(-3.89554531443896569E20+e*(-5.75696351887531875E18-1.238056670725875E16*e))))))))))))))
#define A16   1.26653197415257012E19+e*(2.09399953059891594E21+e*(9.10889795810528434E22+e*(1.63960163245895118E24+e*(1.48019591819210871E25+e*(7.42789224401858187E25+e*(2.19979885688242617E26+e*(3.98058840769200926E26+e*(4.47816195865351041E26+e*(3.1425697955463231E26+e*(1.36178024473674001E26+e*(3.55247020366106089E25+e*(5.32870530549159134E24+e*(4.25081904711579936E23+e*(1.57049964794918696E22+e*(2.0264511586441122E20+3.8379756792502125E17*e)))))))))))))))
#define A17   -4.43286190953399544E20+e*(-8.28945177082857147E22+e*(-4.11156807833097145E24+e*(-8.51681959082844086E25+e*(-8.9426605703698629E26+e*(-5.28429942794582808E27+e*(-1.86982902835006224E28+e*(-4.11362386237013693E28+e*(-5.74697451360533836E28+e*(-5.14202982796267117E28+e*(-2.9383027588358121E28+e*(-1.05685988558916562E28+e*(-2.32509174829616435E27+e*(-2.9808868567899543E26+e*(-2.05578403916548572E25+e*(-6.63156141666285717E23+e*(-7.53586524620779224E21-1.26653197415257012E19*e))))))))))))))))
#define A18   1.64015890652757831E22+e*(3.44433370370791445E24+e*(1.93227120778014001E26+e*(4.56384056694737831E27+e*(5.51464068506141545E28+e*(3.79006214355130008E29+e*(1.5791925598130417E30+e*(4.15102044293713818E30+e*(7.05063031116528617E30+e*(7.83403367907254019E30+e*(5.707653109038565E30+e*(2.70718724539378577E30+e*(8.21180131102781683E29+e*(1.54409939181719633E29+e*(1.71144021260526687E28+e*(1.030544644149408E27+e*(2.92768364815172729E25+e*(2.95228603174964096E23+4.43286190953399544E20*e)))))))))))))))))
#define A19   -6.39661973545755542E23+e*(-1.49894122467555382E26+e*(-9.44332971545598906E27+e*(-2.52271808112895708E29+e*(-3.4757449117776742E30+e*(-2.74899824840597868E31+e*(-1.33220684345828198E32+e*(-4.12349737260896802E32+e*(-8.36827407970643511E32+e*(-1.13045105989016755E33+e*(-1.02278905418634207E33+e*(-6.18524605891345204E32+e*(-2.47409842356538081E32+e*(-6.41432924628061693E31+e*(-1.04272347353330226E31+e*(-1.00908723245158283E30+e*(-5.35122017209172713E28+e*(-1.34904710220799844E27+e*(-1.21535774973693553E25-1.64015890652757831E22*e))))))))))))))))))
#define A20   2.62261409153759772E25+e*(6.81879663799775407E27+e*(4.79361403651242111E29+e*(1.43808421095372633E31+e*(2.24101456206955687E32+e*(2.02098767779363674E33+e*(1.12708928184645126E34+e*(4.05752141464722454E34+e*(9.69628279235549981E34+e*(1.56501406473106313E35+e*(1.72151547120416944E35+e*(1.29283770564739997E35+e*(6.59347229880173987E34+e*(2.25417856369290252E34+e*(5.05246919448409185E33+e*(7.17124659862258199E32+e*(6.11185789655333692E31+e*(2.87616842190745267E30+e*(6.47785680609786637E28+e*(5.24522818307519544E26+6.39661973545755542E23*e)))))))))))))))))))
  double ω = 0;
#if defined (UP_TO_20_TERMS)
  const size_t n_thresholds = 16;
  static const double thresholds[16] = { 10.589, 10.876, 11.22, 11.635, 12.143, 12.771, 13.559, 14.566, 15.884, 17.656, 20.129, 23.743, 29.365, 38.892, 57.148, 99.336 };
  switch (std::upper_bound(thresholds, thresholds + n_thresholds, -h - t + TAU) - thresholds) {
    case 0:  ω = q * (A20 + ω); [[fallthrough]];
    case 1:  ω = q * (A19 + ω); [[fallthrough]];
    case 2:  ω = q * (A18 + ω); [[fallthrough]];
    case 3:  ω = q * (A17 + ω); [[fallthrough]];
    case 4:  ω = q * (A16 + ω); [[fallthrough]];
    case 5:  ω = q * (A15 + ω); [[fallthrough]];
    case 6:  ω = q * (A14 + ω); [[fallthrough]];
    case 7:  ω = q * (A13 + ω); [[fallthrough]];
    case 8:  ω = q * (A12 + ω); [[fallthrough]];
    case 9:  ω = q * (A11 + ω); [[fallthrough]];
    case 10:  ω = q * (A10 + ω); [[fallthrough]];
    case 11:  ω = q * (A9 + ω); [[fallthrough]];
    case 12:  ω = q * (A8 + ω); [[fallthrough]];
    case 13:  ω = q * (A7 + ω); [[fallthrough]];
    case 14:  ω = q * (A6 + ω); [[fallthrough]];
    case 15:  ω = q * (A5 + ω); [[fallthrough]];
    default: ω = A0 + q * (A1 + q * (A2 + q * (A3 + q * (A4 + ω))));
  }
#else
  static const double thresholds[12] = { 12.347, 12.958, 13.729, 14.718, 16.016, 17.769, 20.221, 23.816, 29.419, 38.93, 57.171, 99.347 };
  switch (std::upper_bound(thresholds, thresholds + 12, -h - t + TAU + 0.5) - thresholds) {
    case 0:  ω = q * (A16 + ω); [[fallthrough]];
    case 1:  ω = q * (A15 + ω); [[fallthrough]];
    case 2:  ω = q * (A14 + ω); [[fallthrough]];
    case 3:  ω = q * (A13 + ω); [[fallthrough]];
    case 4:  ω = q * (A12 + ω); [[fallthrough]];
    case 5:  ω = q * (A11 + ω); [[fallthrough]];
    case 6:  ω = q * (A10 + ω); [[fallthrough]];
    case 7:  ω = q * (A9 + ω); [[fallthrough]];
    case 8:  ω = q * (A8 + ω); [[fallthrough]];
    case 9:  ω = q * (A7 + ω); [[fallthrough]];
    case 10:  ω = q * (A6 + ω); [[fallthrough]];
    case 11:  ω = q * (A5 + ω); [[fallthrough]];
    default: ω = A0 + q * (A1 + q * (A2 + q * (A3 + q * (A4 + ω))));
  }
#endif
#endif
  // Note that vega = ∂(b(x,s)/∂s = exp(-(h²+t²)/2)/√(2π)
  const double bx = (t / r) * ω;
  return bx;
}

#if defined( DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS )
double normalised_black_using_erfcx(double h, double t) {
  // Given h = x/s and t = s/2, the normalised Black function can be written as
  //
  //     b(x,s,θ)  =  θ · [ Φ(θ·(x/s+s/2))·exp(x/2)  -   Φ(θ·(x/s-s/2))·exp(-x/2) ]
  //               =  θ · [ Φ(h+θ·t)·exp(x/2)      -   Φ(h-θ·t)·exp(-x/2) ]
  //               =  Φ(h+t)·exp(θ·x/2)    -   Φ(h-t)·exp(-θ·x/2)
  //               =  Φ(h+t)·exp(h·t)      -   Φ(h-t)·exp(-h·t)                       (*)
  // with
  //                h  =  θ·x/s   and   t  =  s/2 .
  // 
  // ⟹  θ can be subsumed into x (or h)
  //
  // It is mentioned in section 4 (and discussion of figures 2 and 3) of George Marsaglia's article "Evaluating the
  // Normal Distribution" (available at http://www.jstatsoft.org/v11/a05/paper) that the error of any cumulative normal
  // function Φ(z) is dominated by the hardware (or compiler implementation) accuracy of exp(-z²/2) which is not
  // reliably more than 14 digits when z is large. The accuracy of Φ(z) typically starts coming down to 14 digits when
  // z is around -8. For the (normalised) Black function, as above in (*), this means that we are subtracting two terms
  // that are each products of terms with about 14 digits of accuracy. The net result, in each of the products, is even
  // less accuracy, and then we are taking the difference of these terms, resulting in even less accuracy. When we are
  // using the asymptotic expansion asymptotic_expansion_of_scaled_normalised_black() invoked in the second branch at the
  // beginning of this function, we are using only *one* exponential instead of 4, and this improves accuracy. It
  // actually improves it a bit more than you would expect from the above logic, namely, almost the full two missing
  // digits (in 64 bit IEEE floating point).  Unfortunately, going higher order in the asymptotic expansion will not
  // enable us to gain more accuracy (by extending the range in which we could use the expansion) since the asymptotic
  // expansion, being a divergent series, can never gain 16 digits of accuracy for z=-8 or just below. The best you can
  // get is about 15 digits (just), for about 35 terms in the series (26.2.12), which would result in an prohibitively
  // long expression in function asymptotic expansion asymptotic_expansion_of_scaled_normalised_black(). In this last branch,
  // here, we therefore take a different tack as follows.
  //     The "scaled complementary error function" is defined as erfcx(z) = exp(z²)·erfc(z). Cody's implementation of this
  // function as published in "Rational Chebyshev approximations for the error function", W. J. Cody, Math. Comp., 1969, pp.
  // 631-638, uses rational functions that theoretically approximates erfcx(x) to at least 18 significant decimal digits,
  // *without* the use of the exponential function when x>4, which translates to about z<-5.66 in Φ(z). To make use of it,
  // we write
  //             Φ(z) = exp(-z²/2)·erfcx(-z/√2)/2
  //
  // to transform the normalised black function to
  //
  //   b   =  half · exp(-half(h²+t²)) · [ erfcx(-(h+t)/√2) -  erfcx(-(h-t)/√2) ]
  //
  // which now involves only one exponential, instead of three, when |h|+|t| > 5.66 , and the difference inside the
  // square bracket is between the evaluation of two rational functions, which, typically, according to Marsaglia,
  // retains the full 16 digits of accuracy (or just a little less than that).
  //
  const double b = 0.5 * exp(-0.5 * (h * h + t * t)) * (erfcx_cody(-(1 / SQRT_TWO) * (h + t)) - erfcx_cody(-(1 / SQRT_TWO) * (h - t)));
  return fabs(std::max(b, 0.0));
}
#endif

inline double Yprime_tail_expansion_rational_function_part(double w) {
  return w * (-2.9999999999994663866 + w * (-1.7556263323542206288E2 + w * (-3.4735035445495633334E3 + w * (-2.7805745693864308643E4 + w * (-8.3836021460741980839E4 - 6.6818249032616849037E4 * w))))) / (1 + w * (6.3520877744831739102E1 + w * (1.4404389037604337538E3 + w * (1.4562545638507033944E4 + w * (6.6886794165651675684E4 + w * (1.2569970380923908488E5 + 6.9286518679803751694E4 * w))))));
}

//
// Y'(h) = 1+h·Y(h) avoiding subtractive cancellation.
//
double Yprime(double h) {
  // We copied the thresholds of -0.46875 and -4 from Cody.
  if (h < -4) {
    // Nonlinear-Remez optimized minimax rational function of order (5,6) for g(w) := (Y'(h)/h²-1)/h² with w:=1/h².
    // The relative accuracy of Y'(h) ≈ w·(1+w·g(w)) is better than 9.8E-17 (in perfect arithmetic) on h in [-∞,-4] (i.e., on w in [0,1/16]).
    const double w = 1 / (h * h);
    return w * (1 + Yprime_tail_expansion_rational_function_part(w));
  }
  if (h <= -0.46875)
    // Remez-optimized minimax rational function of order (7,7) of relative accuracy better than 1.6E-16 (in perfect arithmetic) on h in [-4,-0.46875].
    return (1.0000000000594317229 - h * (6.1911449879694112749E-1 - h * (2.2180844736576013957E-1 - h * (4.5650900351352987865E-2 - h * (5.545521007735379052E-3 - h * (3.0717392274913902347E-4 - h * (4.2766597835908713583E-8 + 8.4592436406580605619E-10 * h))))))) / (1 - h * (1.8724286369589162071 - h * (1.5685497236077651429 - h * (7.6576489836589035112E-1 - h * (2.3677701403094640361E-1 - h * (4.6762548903194957675E-2 - h * (5.5290453576936595892E-3 - 3.0822020417927147113E-4 * h)))))));
  return 1 + h * SQRT_PI_OVER_TWO * erfcx_cody(-(1 / SQRT_TWO) * h);
}

// Calculation of the 'scaled normalised Black' function (in analogy to the 'scaled' complementary error function)
// 
//             bx :=  b  /  ∂(b(x,s)/∂s
// 
// with h := θ·x/s, t := s/2, and
//
//              b  =  Φ(h+t)·exp(h·t) - Φ(h-t)·exp(-h·t)
//
//                    exp(-(h²+t²)/2)
//                 =  --------------- ·  [ Y(h+t) - Y(h-t) ]
//                        √(2π)
// and
//        ∂(b(x,s)/∂s = exp(-(h²+t²)/2)/√(2π)
// and
//           Y(z) := Φ(z)/φ(z)
// whence
//             bx  =  Y(h+t) - Y(h-t)
//
// using an expansion for small t to twelfth order in t.
// Theoretically accurate to (better than) precision  ε = 2.23E-16  when  h<=0  and  t < τ  with  τ := 2·ε^(1/16) ≈ 0.21.
// The main bottleneck for precision is the coefficient a:=1+h·Y(h) when |h|>1 .
//
// Note that vega = ∂(b(x,s)/∂s = exp(-(h²+t²)/2)/√(2π)
//
double small_t_expansion_of_scaled_normalised_black(double h, double t) {
  // Y(h) := Φ(h)/φ(h) = √(π/2)·erfcx(-h/√2)
  // a := 1+h·Y(h)  --- Note that due to h<0, and h·Y(h) -> -1 (from above) as h -> -∞, we also have that a>0 and a -> 0 as h -> -∞
  // g := h² , ζ := t²
#if defined( USE_ORIGINAL_TAYLOR_EXPANSION )
  const double a = 1 + h * SQRT_PI_OVER_TWO * erfcx_cody(-(1 / SQRT_TWO) * h), g = h * h, ζ = t * t;
  return 2 * t * (a + ζ * ((-1 + 3 * a + a * g) / 6 + ζ * ((-7 + 15 * a + g * (-1 + 10 * a + a * g)) / 120 + ζ * ((-57 + 105 * a + g * (-18 + 105 * a + g * (-1 + 21 * a + a * g))) / 5040 + ζ * ((-561 + 945 * a + g * (-285 + 1260 * a + g * (-33 + 378 * a + g * (-1 + 36 * a + a * g)))) / 362880 + ζ * ((-6555 + 10395 * a + g * (-4680 + 17325 * a + g * (-840 + 6930 * a + g * (-52 + 990 * a + g * (-1 + 55 * a + a * g))))) / 39916800 + ((-89055 + 135135 * a + g * (-82845 + 270270 * a + g * (-20370 + 135135 * a + g * (-1926 + 25740 * a + g * (-75 + 2145 * a + g * (-1 + 78 * a + a * g)))))) * ζ) / 6227020800.0))))));
#else
  const double a = Yprime(h), h² = h * h, t² = t * t;
# define B0 2*a
# define B1 (-1+a*(3+h²))/3
# define B2 (-7-h²+a*(15+h²*(10+h²)))/60
# define B3 (-57+(-18-h²)*h²+a*(105+h²*(105+h²*(21+h²))))/2520
# define B4 (-561+h²*(-285+(-33-h²)*h²)+a*(945+h²*(1260+h²*(378+h²*(36+h²)))))/181440
# define B5 (-6555+h²*(-4680+h²*(-840+(-52-h²)*h²))+a*(10395+h²*(17325+h²*(6930+h²*(990+h²*(55+h²))))))/19958400
# define B6 (-89055+h²*(-82845+h²*(-20370+h²*(-1926+(-75-h²)*h²)))+a*(135135+h²*(270270+h²*(135135+h²*(25740+h²*(2145+h²*(78+h²)))))))/3113510400
  return t * (B0 + t² * (B1 + t² * (B2 + t² * (B3 + t² * (B4 + t² * (B5 + B6 * t²))))));
#endif
}

#if defined( DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS )
double normalised_black_using_norm_cdf(double θx, double s) {
  //     b(x,s,θ)  =  θ · [ Φ(θ·(x/s+s/2))·exp(x/2)  -   Φ(θ·(x/s-s/2))·exp(-x/2) ]
  //               =  θ · [ Φ(h+θ·t)·exp(x/2)      -   Φ(h-θ·t)·exp(-x/2) ]
  //               =  Φ(h+t)·exp(θ·x/2)    -   Φ(h-t)·exp(-θ·x/2)
  // with
  //              h  =  θ·x/s   and   t  =  s/2
  const double h = θx / s, t = 0.5 * s, b_max = exp(0.5 * θx), b = norm_cdf(h + t) * b_max - norm_cdf(h - t) / b_max;
  return fabs(std::max(b, 0.0));
}
#endif

//
// Introduced on 2017-02-18
//
//     b(x,s,θ)  =  θ · [ Φ(θ·(x/s+s/2))·exp(x/2)  -   Φ(θ·(x/s-s/2))·exp(-x/2) ]
//               =  θ · [ Φ(h+θ·t)·exp(x/2)      -   Φ(h-θ·t)·exp(-x/2) ]
//               =  Φ(h+t)·exp(θ·x/2)    -   Φ(h-t)·exp(-θ·x/2)
// with
//              h  =  θ·x/s ,       t  =  s/2 ,
//
// ⟹  θ can be subsumed into x (or h)
//
//     b(x,s)  =  half · exp(-u²-v²) · [ erfcx(u+v) -  erfcx(u-v) ]
//             =  half · [ exp(θx/2)·erfc(u+v)     -  exp(-θx/2)·erfc(u-v)    ]
//             =  half · [ exp(θx/2)·erfc(u+v)     -  exp(-u²-v²)·erfcx(u-v) ]
//             =  half · [ exp(-u²-v²)·erfcx(u+v) -  exp(-θx/2)·erfc(u-v)    ]
// and
//              u  = -h/√2  and   v  =  -t/√2 .
//
// Cody's erfc() and erfcx() functions each, for some values of their argument, involve the evaluation
// of the exponential function exp(). The normalised Black function requires additional evaluation(s)
// of the exponential function irrespective of which of the above formulations is used. However, the total
// number of exponential function evaluations can be minimised by a judicious choice of one of the above
// formulations depending on the input values and the branch logic in Cody's erfc() and erfcx().
//
double normalised_black_with_optimal_use_of_codys_functions(double θx, double s) {
  const double codys_threshold = 0.46875, h = θx / s, t = 0.5 * s, q₁ = -(1 / SQRT_TWO) * (h + t), q₂ = -(1 / SQRT_TWO) * (h - t);
  double two_b;
  // Note: the absence of fabs() around q₁ and q₂ is no accident.
  if (q₁ < codys_threshold)
    if (q₂ < codys_threshold)
      two_b = exp(0.5 * θx) * erfc_cody(q₁) - exp(-0.5 * θx) * erfc_cody(q₂);
    else
      two_b = exp(0.5 * θx) * erfc_cody(q₁) - exp(-0.5 * (h * h + t * t)) * erfcx_cody(q₂);
  else
    if (q₂ < codys_threshold)
      two_b = exp(-0.5 * (h * h + t * t)) * erfcx_cody(q₁) - exp(-0.5 * θx) * erfc_cody(q₂);
    else
      two_b = exp(-0.5 * (h * h + t * t)) * (erfcx_cody(q₁) - erfcx_cody(q₂));
  return std::max(0.5 * two_b, 0.0);
}

// ∂b(x,s)/∂s = b'(s) = exp(-half·((x/s)²+(s/2)²) / √(2π)
inline double normalised_vega(double x, double s) {
#if defined( SAFEGUARD_VEGA_AGAINST_ZERO_VOLATILITY )
  const double ax = fabs(x);
  if (ax <= 0)
    return (1 / SQRT_TWO_PI) * exp(-0.125 * s * s);
  if (s <= 0 || s <= ax * SQRT_DBL_MIN) {
#if defined(LOG_VEGA_OVERFLOW)
    printf("normalised_vega(x = %g, s = %g)\n", x, s);
#endif
    return 0;
  }
#else
  assert(s > 0);
#endif
  const double h = x / s, t = 0.5 * s;
  return (1 / SQRT_TWO_PI) * exp(-0.5 * (h * h + t * t));
}

inline double inv_normalised_vega(double x, double s) {
#if defined( SAFEGUARD_VEGA_AGAINST_ZERO_VOLATILITY )
  const double ax = fabs(x);
  if (ax <= 0)
    return SQRT_TWO_PI * exp(0.125 * s * s);
  if (s <= 0 || s <= ax * SQRT_DBL_MIN) {
#if defined(LOG_VEGA_OVERFLOW)
    printf("inv_normalised_vega(x = %g, s = %g)\n", x, s);
#endif
    return DBL_MAX;
  }
#else
  assert(s > 0);
#endif
  const double h = x / s, t = 0.5 * s;
  return SQRT_TWO_PI * exp(0.5 * (h * h + t * t));
}

inline double ln_normalised_vega(double x, double s) {
#if defined( SAFEGUARD_VEGA_AGAINST_ZERO_VOLATILITY )
  if (fabs(x) <= 0)
    return (-(LN_TWO_PI / 2) - 0.125 * (s * s));
  if (s <= 0) {
#if defined(LOG_VEGA_OVERFLOW)
    printf("ln_normalised_vega(x = %g, s = %g)\n", x, s);
#endif
    return -DBL_MAX;
  }
#else
  assert(s > 0);
#endif
  const double h = x / s, t = 0.5 * s;
  return -(LN_TWO_PI * 0.5) - 0.5 * (h * h + t * t);
}

#define SPECIALISE_EXACT_ATM_CASE

//  b(x,s,θ)  =  θ·[ exp(θx/2)·Φ(θ·(x/s+s/2)) - exp(-θx/2)·Φ(θ·(x/s-s/2)) ]
//            =  exp(θx/2)·Φ(θx/s+s/2) - exp(-θx/2)·Φ(θx/s-s/2)
//  ⟹  θ can be subsumed into x.
// NOTE: this function requires s > 0 and θx ≤ 0 or even θx < 0 when SPECIALISE_EXACT_ATM_CASE is defined.
double normalised_black(double θx, double s) {
#if defined( SPECIALISE_EXACT_ATM_CASE )
  assert(θx < 0);
#else
  assert(θx <= 0);
#endif
//  if (VOLATILITY_IS_CONSIDERED_ZERO(θx, s)) return 0;
  assert(s > 0);
  if (IS_REGION_I(θx, s))
    return asymptotic_expansion_of_scaled_normalised_black(θx / s, 0.5 * s) * normalised_vega(θx, s);
  if (IS_REGION_II(θx, s))
    return small_t_expansion_of_scaled_normalised_black(θx / s, 0.5 * s) * normalised_vega(θx, s);
#if defined( DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS )
  if (IS_REGION_III(θx, s))
    return normalised_black_using_norm_cdf(θx, s);
  return normalised_black_using_erfcx(θx / s, 0.5 * s);
#else
  return normalised_black_with_optimal_use_of_codys_functions(θx, s);
#endif
}

// NOTE: this function requires s > 0 and θx ≤ 0 or even θx < 0 when SPECIALISE_EXACT_ATM_CASE is defined.
std::tuple<double, double> scaled_normalised_black_and_ln_vega(double θx, double s) {
#if defined( SPECIALISE_EXACT_ATM_CASE )
  assert(θx < 0);
#else
  assert(θx <= 0);
#endif
//  if (VOLATILITY_IS_CONSIDERED_ZERO(θx, s)) return { 0, -DBL_MAX };
  assert(s > 0);
  if (IS_REGION_I(θx, s))
    return { asymptotic_expansion_of_scaled_normalised_black(θx / s, 0.5 * s), ln_normalised_vega(θx, s) };
  if (IS_REGION_II(θx, s))
    return { small_t_expansion_of_scaled_normalised_black(θx / s, 0.5 * s), ln_normalised_vega(θx, s) };
  const double ln_vega = ln_normalised_vega(θx, s);
#if defined( DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS )
  if (IS_REGION_III(θx, s))
    return { normalised_black_using_norm_cdf(θx, s) * exp(-ln_vega), ln_vega };
  return { normalised_black_using_erfcx(θx / s, 0.5 * s) * exp(-ln_vega), ln_vega };
#else
  return { normalised_black_with_optimal_use_of_codys_functions(θx, s) * exp(-ln_vega), ln_vega };
#endif
}

// ∂²b(x,s)/∂s²
double normalised_volga(double x, double s) {
  const double ax = fabs(x);
  if (ax <= 0) return (1 / SQRT_TWO_PI) * exp(-0.125 * s * s);
  if (s <= 0 || s <= ax * SQRT_DBL_MIN) return 0;
  const double h = x / s, t = 0.5 * s, h² = h * h, t² = t * t;
  return (1 / SQRT_TWO_PI) * exp(-0.5 * (h² + t²)) * (h² - t²) / s;
}

#ifdef COMPUTE_LOWER_MAP_DERIVATIVES_INDIVIDUALLY
double f_lower_map(const double x, const double s) {
#if defined( POSITIVE_DENORMALISATION_CUTOFF )
  if (is_below_horizon(x)) return 0;
  if (is_below_horizon(s)) return 0;
#endif
  const double z = SQRT_ONE_OVER_THREE * fabs(x) / s, Phi = norm_cdf(-z);
  return TWO_PI_OVER_SQRT_TWENTY_SEVEN * fabs(x) * (Phi * Phi * Phi);
}
double d_f_lower_map_d_beta(const double x, const double s) {
#if defined( POSITIVE_DENORMALISATION_CUTOFF )
  if (is_below_horizon(s)) return 1;
#endif
  const double z = SQRT_ONE_OVER_THREE * fabs(x) / s, y = z * z, Phi = norm_cdf(-z);
  return TWO_PI * y * (Phi * Phi) * exp(y + 0.125 * s * s);
}
double d2_f_lower_map_d_beta2(const double x, const double s) {
  const double ax = fabs(x), z = SQRT_ONE_OVER_THREE * ax / s, y = z * z, s2 = s * s, Phi = norm_cdf(-z), phi = norm_pdf(z);
  return PI_OVER_SIX * y / (s2 * s) * Phi * (8 * SQRT_THREE * s * ax + (3 * s2 * (s2 - 8) - 8 * x * x) * Phi / phi) * exp(2 * y + 0.25 * s2);
}
void compute_f_lower_map_and_first_two_derivatives(const double x, const double s, double& f, double& fp, double& fpp) {
  f = f_lower_map(x, s);
  fp = d_f_lower_map_d_beta(x, s);
  fpp = d2_f_lower_map_d_beta2(x, s);
}
#else
inline void compute_f_lower_map_and_first_two_derivatives(const double x, const double s, double& f, double& fp, double& fpp) {
  const double ax = fabs(x), z = SQRT_ONE_OVER_THREE * ax / s, y = z * z, s2 = s * s, Phi = 0.5 * erfc_cody((1 / SQRT_TWO) * z) /* = norm_cdf(-z) */, phi = norm_pdf(z);
  fpp = PI_OVER_SIX * y / (s2 * s) * Phi * (8 * SQRT_THREE * s * ax + (3 * s2 * (s2 - 8) - 8 * x * x) * Phi / phi) * exp(2 * y + 0.25 * s2);
#if defined( POSITIVE_DENORMALISATION_CUTOFF )
  if (is_below_horizon(s)) {
    fp = 1;
    f = 0;
  } else
#endif
  {
    const double Phi2 = Phi * Phi;
    fp = TWO_PI * y * Phi2 * exp(y + 0.125 * s * s);
#if defined( POSITIVE_DENORMALISATION_CUTOFF )
    if (is_below_horizon(x)) f = 0; else
#endif
      f = TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * (Phi2 * Phi);
  }
}
#endif

// Formula (4.38)
inline double inverse_f_lower_map(const double x, const double f) {
  // Caution: this can result in unnecessary underflow when f ≈ DBL_MIN and |x| is large. It also triggers 'unsafe-math-optimizations' issues with g++, exposed when -Ofast or -fast-math is used.
  //   return is_below_horizon(f) ? 0 : fabs(x / (SQRT_THREE * inverse_norm_cdf(std::pow(f / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * fabs(x)), 1. / 3.))));
  // Second caution: do *NOT* attempt to contract the two cubic roots into one. This leads to catastrophic underflow when f ≈ DBL_MIN and |x| > 1.
#if defined( POSITIVE_DENORMALISATION_CUTOFF )
  if (is_below_horizon(f)) return DBL_MIN;
#endif
  return fabs(x / (SQRT_THREE * inverse_norm_cdf(SQRT_THREE_OVER_THIRD_ROOT_TWO_PI * std::cbrt(f) / std::cbrt(fabs(x)))));
}

#ifdef COMPUTE_UPPER_MAP_DERIVATIVES_INDIVIDUALLY
double f_upper_map(const double s) {
  return norm_cdf(-0.5 * s);
}
double d_f_upper_map_d_beta(const double x, const double s) {
#if defined( POSITIVE_DENORMALISATION_CUTOFF )
  if (is_below_horizon(x)) return -0.5;
#endif
  return -0.5 * exp(0.5 * square(x / s));
}
double d2_f_upper_map_d_beta2(const double x, const double s) {
#if defined( POSITIVE_DENORMALISATION_CUTOFF )
  if (is_below_horizon(x)) return 0;
#endif
  const double w = square(x / s);
  return SQRT_PI_OVER_TWO * exp(w + 0.125 * s * s) * w / s;
}
void compute_f_upper_map_and_first_two_derivatives(const double x, const double s, double& f, double& fp, double& fpp) {
  f = f_upper_map(s);
  fp = d_f_upper_map_d_beta(x, s);
  fpp = d2_f_upper_map_d_beta2(x, s);
}
#else
inline void compute_f_upper_map_and_first_two_derivatives(const double x, const double s, double& f, double& fp, double& fpp) {
  f = 0.5 * erfc_cody((0.5 / SQRT_TWO) * s) /* = norm_cdf(-0.5 * s) */;
#if defined( SAFEGUARD_AGAINST_NON_POSITIVE_VOLATILITY )
  if (is_below_horizon(x)) { fp = -0.5; fpp = 0; return; }
#endif
  const double w = square(x / s);
  fp = -0.5 * exp(0.5 * w);
  fpp = SQRT_PI_OVER_TWO * exp(w + 0.125 * s * s) * w / s;
}
#endif

inline double inverse_f_upper_map(double f) { return -2. * inverse_norm_cdf(f); }

// Extensive testing indicates that this is never needed.
#if defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
inline double take_step(const double x_min, const double x_max, double x, double& dx) {
  const double new_x = x + dx;
  if ( new_x > x_max ) {
    dx = x_max-x;
    return x_max;
  }
  if ( new_x < x_min ) {
    dx = x_min-x;
    return x_min;
  }
  return new_x;
}
#endif

//
// Introduced on 2023-12-10
// 
//    b̄(x,s,θ)          :=   bₘₐₓ(x,θ)   -  b(x,s,θ)
//                       =   exp(θ·x/2)  -  θ·[ exp(θx/2)·Φ(θ·(x/s+s/2)) - exp(-θx/2)·Φ(θ·(x/s-s/2)) ]            using  bₘₐₓ = exp(θ·x/2)
//
//                           ⎧   exp(θx/2)·[1-Φ(x/s+s/2)] + exp(-θx/2)·Φ(x/s-s/2)                                  when θ = +1
//                       =   ⎨ 
//                           ⎩   exp(θx/2)·Φ(-x/s-s/2)    + exp(-θx/2)·[1-Φ(-x/s+s/2)]                             when θ = -1
// using  1-Φ(z) = Φ(-z)
//                       =   exp(θx/2)·Φ(-x/s-s/2) + exp(-θx/2)·Φ(x/s-s/2)                            (1)          for both θ = ±1
// 
// Note: no subtractive cancellation, and no dependency on θ = ±1 !
// 
//    b̄(x,s)    =   half · [ erfc((s/2+x/s)/√2)·exp(θx/2) + erfc((s/2-x/s)/√2)·exp(-θx/2) ]              (2)          using  Φ(z) = erfc(-z/√2)/2
//
// using erfc(z) = erfcx(z)·exp(-z²)
// 
//    b̄(x,s)    =   half · ( erfcx((s/2+x/s)/√2)·exp(-half((x/s)²+x+(s/2)²))·exp(θx/2) + erfcx((s/2-x/s)/√2)·exp(-half((x/s)²-x+(s/2)²))·exp(-θx/2) )
//
//              =   half · ( erfcx((s/2+x/s)/√2) + erfcx((s/2-x/s)/√2) ) · exp(-half((x/s)²+(s/2)²))
//
//              =   √(π/2) · [ erfcx((s/2+x/s)/√2) + erfcx((s/2-x/s)/√2) ] · ∂b(x,s)/∂s
//
inline double complementary_normalised_black(double h /* = x/s */, double t /* = s/2 */) {
  //    b̄(x,s)    =   half · ( erfcx((s/2+x/s)/√2) + erfcx((s/2-x/s)/√2) ) · exp(-half((x/s)²+(s/2)²))
  //    b̄(x,s)    =   half · ( erfcx((t+h)/√2) + erfcx((t-h)/√2) ) · exp(-half(t²+h²))
  return 0.5 * (erfcx_cody((t + h) * (1 / SQRT_TWO)) + erfcx_cody((t - h) * (1 / SQRT_TWO))) * exp(-0.5 * (t * t + h * h));
}

// f(x) := 1 - erfcx(x)  ≈  2/√π·x - x² + 4/(3√π)·x³ + ... for small x.
inline double one_minus_erfcx(double x) {
  if (x < -1. / 5. || x > 1. / 3.)
    return 1 - erfcx_cody(x);
  // Remez-optimized minimax rational function of order (4,5) for g(x) := (2/√π-f(x)/x)/x.   2/√π = 1.128379167095512573896.
  // The relative accuracy of f(x) ≈ x·(2/√π-x·g(x)) is better than 2.5E-17 (in perfect arithmetic) on x in [-1/5,1/3].
  return x * (1.128379167095512573896 - x * (1.0000000000000002 + x * (1.1514967181784756 + x * (5.7689001208873741E-1 + x * (1.4069188744609651E-1 + 1.4069285713634565E-2 * x)))) / (1 + x * (1.9037494962421563 + x * (1.5089908593742723 + x * (6.2486081658640257E-1 + x * (1.358008134514386E-1 + 1.2463320728346347E-2 * x))))));
}

#define USE_UNIVARIATE_RATIONAL_FUNCTION_APPROXIMATIONS_FOR_BL_AND_BU_OVER_BMAX

#if defined( USE_UNIVARIATE_RATIONAL_FUNCTION_APPROXIMATIONS_FOR_BL_AND_BU_OVER_BMAX )

//  s_c        :=   √(2|x|)
//  b_c(x)     :=   b(-|x|,s_c)
//  vega(x,s)  :=   ∂b(x,s)/∂s
//  sₗ          :=   s_c - b_c(x)/vega(x,s_c)
//              =   √(-2x) - √(π/2) · (1 - erfcx(√(-x)))
//  bₗ          :=   b(-|x|,sₗ)
//  sᵤ         :=  s_c + (bₘₐₓ-b_c(x))/vega(x,s_c)
//              =  √(-2x) + √(π/2) · (1 + erfcx(√(-x)))
//  bᵤ         :=  b(-|x|,sᵤ)
//  bₘₐₓ       :=  exp(-|x|/2)
// 
// Note that bₗ=bₗ(x), bᵤ=bᵤ(x), and bₘₐₓ=bₘₐₓ(x) are all univariate functions.

// Within the "Let's Be Rational" algorithm, the functions bₗ_over_bₘₐₓ and bᵤ_over_bₘₐₓ are only
// used as an acceleration for bₗ and bᵤ which in turn are only used to determine branches and for
// the initial guess interpolation logic. Thus, if their accuracy is only (relative) 1.E-15 or even
// 1.E-14 is perfectly sufficient. Recall that the initial guess accuracy can be as low as 1E-4, so,
// any amount of noise that is 1E-10 smaller than that is of zero subsequent significance.

// Nevertheless, we prefer the four-branch versions further down since they are slightly faster than the
// two-branch versions (on the hardware used for testing, i.e., a 12th Gen Intel(R) Core(TM) i5-12500H).

#if defined( BL_AND_BU_OVER_BMAX_IN_TWO_BRANCHES )

inline double bₗ_over_bₘₐₓ(const double s_c /* = √|2x| */) {
  if (s_c >= 2.82842712474619 /* = √8 */) {
    // x = -1/(2·y²), y = 1/√|2x|, output is f(x) = bₗ(x)/bₘₐₓ(x).
    // Nonlinear-Remez optimized minimax rational function of order (7,7) for g(y) := ((f - Φ(-√(π/2)) )·4·exp(π/4)/y+√(π/2))/y .
    // Accuracy better than 2.3E-16 in perfect arithmetic on x in [-∞,-4], i.e., y in [0,1/√8].
    const double y = 1 / s_c, g = (4.6053948312498956145E-2 + y * (2.2003742013556747214 + y * (1.7600379430086971755E1 + y * (9.1467986769463622162E1 + y * (2.9140932094514209454E2 + y * (5.7000507935678940271E2 + y * (5.6859094292134607288E2 - 2.2892636835707606708E0 * y))))))) / (1 + y * (8.5913225035477097129 + y * (5.0371467968780246369E1 + y * (1.8369076530154652307E2 + y * (4.6283633310610371326E2 + y * (7.8579932950540139912E2 + y * (8.0484923741700038479E2 + 4.4439123028067025538E2 * y)))))));
    // f = Φ(-√(π/2)) + exp(-π//4)/4 · y · ( -√(π/2) + y · g )
    return 0.1050457027219686376856 + 0.1139845319414990591915 * y * (-1.253314137315500251208 + y * g);
  }
  // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x).
  // For small |x|, i.e., small y, we have f ≈ (exp(-1/π)/4-Φ(-√(π/2))/2)·y² + exp(-1/π)/(3·√(2π))·y³ + O(y⁴).
  //   c₂ := (exp(-1/π)/4-Φ(-√(π/2))/2) =  0.07560996640296361767172
  //   c₃ := exp(-1/π)/(3·√(2π))        = -0.09672719281339436290858
  // Nonlinear-Remez optimized minimax rational function of order (7,7) for g(y) := ((f/y²-c₂)/y-c₃)/y .
  // Accuracy better than 1E-15 in perfect arithmetic on x in [-4,0], i.e., y in [0,√8].
  const double g = (8.074107237288007116E-2 + s_c * (2.1690496514582263229E-1 + s_c * (2.3522503078524868067E-1 + s_c * (1.3517614695211948957E-1 + s_c * (4.4661233146407726342E-2 + s_c * (8.2133824762408676406E-3 + s_c * (6.7102169991726196873E-4 - 1.7070712794798401108E-10 * s_c))))))) / (1 + s_c * (3.331190548448924999 + s_c * (4.7356539912342037297 + s_c * (3.7552272118554961136 + s_c * (1.8061741297788241057 + s_c * (5.3232597992823920757E-1 + s_c * (9.0336639624394764651E-2 + 6.9371970943430933891E-3 * s_c)))))));
  // f = y²·(c₂+y·(c₃+y·g))
  return (s_c * s_c) * (0.07560996640296361767172 + s_c * (-0.09672719281339436290858 + s_c * g));
}

inline double bᵤ_over_bₘₐₓ(const double s_c /* = √|2x| */) {
  if (s_c >= 2.449489742783178098197 /* = √6 */) {
    // x = -1/(2·y²), y = 1/√|2x|, output is f(x) = bᵤ(x)/bₘₐₓ(x).
    // Nonlinear-Remez optimized minimax rational function of order (7,7) for g(y) := ((f - Φ(√(π/2)) )·4·exp(π/4)/y+√(π/2))/y .
    // Accuracy better than 7.8E-17 in perfect arithmetic on x in [-∞,-3], i.e., y in [0,1/√6].
    const double y = 1 / s_c, g = (-4.6053948172126087243E-2 + y * (1.3751630820772591756 + y * (1.455319886249397753E1 + y * (8.4880892200802385788E1 + y * (2.9836801628056627321E2 + y * (6.1696928351291694236E2 + y * (6.5892809576774076373E2 - 1.2291897122716543832E0 * y))))))) / (1 + y * (9.3270349037904060715 + y * (5.4363781465880728853E1 + y * (2.0304204599521774403E2 + y * (5.0796471791232277633E2 + y * (8.6988303136901841628E2 + y * (8.8812383339606771272E2 + 5.2060847522792563067E2 * y)))))));
    // f = Φ(√(π/2)) + exp(-π//4)/4 · y · ( -√(π/2) + y · g )
    return 0.8949542972780313623144 + 0.1139845319414990591915 * y * (-1.253314137315500251208 + y * g);
  }
  // x = -y²/2, y = √|2x| = s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x).
  // For small |x|, i.e., small y, we have f ≈ 1-2·Φ(-√(π/2)) + (exp(-π/4)/4-Φ(-√(π/2))/2)·y² + O(y³).
  //   c₀ := 1-2·Φ(-√(π/2))             = 0.7899085945560627246288
  //   c₂ := (exp(-π/4)/4-Φ(-√(π/2))/2) = 0.06146168058051474034868
  // Remez-optimized minimax rational function of order (7,7) for g(y) := ((f - c₀)/y²-c₂)/y .
  // Accuracy better than 2.1E-17 in perfect arithmetic on x in [-3,0], i.e., y in [0,√6].
  const double g = (-6.0630998803348508652E-2 + s_c * (-1.3448643785893710638E-1 + s_c * (-1.4163681164247209788E-1 + s_c * (-8.9763830861375452184E-2 + s_c * (-3.5121337410416895427E-2 + s_c * (-8.1438128785484908289E-3 + s_c * (-8.7339910261568864981E-4 - 3.3867568170011767283E-9 * s_c))))))) / (1 + s_c * (2.7220033406555055453 + s_c * (3.6503350360158844006 + s_c * (3.0186389537663895754 + s_c * (1.6527347941968487339 + s_c * (5.9591616493512210562E-1 + s_c * (1.3248016238920728894E-1 + 1.4212067435291777677E-2 * s_c)))))));
  // f = c₀ + y²·(c₂+y·g)
  return 0.7899085945560627246288 + (s_c * s_c) * (0.06146168058051474034868 + s_c * g);
}

#else

inline double bₗ_over_bₘₐₓ(const double s_c /* = √|2x| */) {
  // Four branches.
  //                    │     I.          │         II.           │        III.         │        IV.     │
  //                    │    ←-→          │         ←-→           │        ←-→          │        ←-→     │
  //  In |x|:           [  0  ,  0.252               ,  3.45                ,  27                 ,  ∞   )
  //  In s_c = √|2x| :  [  0  ,  0.7099295739719539  ,  2.6267851073127395  ,  7.348469228349534  ,  ∞   )
#if !defined(BL_BRANCHES_UNORDERED)
  if (s_c < 2.6267851073127395) {
    if (s_c < 0.7099295739719539) {
      // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x).
      // For small |x|, i.e., small y, we have f ≈ (exp(-1/π)/4-Φ(-√(π/2))/2)·y² + exp(-1/π)/(3·√(2π))·y³ + O(y⁴).
      //   c₂ := (exp(-1/π)/4-Φ(-√(π/2))/2) =  0.07560996640296361767172
      //   c₃ := exp(-1/π)/(3·√(2π))        = -0.09672719281339436290858
      // Nonlinear-Remez optimized minimax rational function of order (5,4) for g(y) := ((f/y²-c₂)/y-c₃)/y .   f = y²·(c₂+y·(c₃+y·g)) .
      const double g = (8.0741072372882856924E-2 + s_c * (9.8078911786358897272E-2 + s_c * (3.9760631445677058375E-2 + s_c * (5.9716928459589189876E-3 + s_c * (-6.4036399341479799981E-6 + 4.5425102093616062245E-7 * s_c))))) / (1 + s_c * (1.8594977672287664353 + s_c * (1.3658801475711790419 + s_c * (4.6132707108655653215E-1 + 6.1254597049831720643E-2 * s_c))));
      // Branch I. Accuracy better than 7.43E-17 in perfect arithmetic.
      return (s_c * s_c) * (0.07560996640296361767172 + s_c * (s_c * g - 0.09672719281339436290858));
    }
    // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,6) for g(y) := bᵤ(-y²/2)/bₘₐₓ(-y²/2).
    // Branch II. Accuracy better than 8.77E-17 in perfect arithmetic.
    return (1.9795737927598581235E-9 + s_c * (-2.7081288564685588037E-8 + s_c * (7.5610142272549044609E-2 + s_c * (6.917130174466834016E-2 + s_c * (2.9537058950963019803E-2 + s_c * (6.5849252702302307774E-3 + 6.9711400639834715731E-4 * s_c)))))) / (1 + s_c * (2.1941448525586579756 + s_c * (2.1297103549995181357 + s_c * (1.1571483187179784072 + s_c * (3.7831622253060456794E-1 + s_c * (7.1714862448829349869E-2 + 6.6361975827861200167E-3 * s_c))))));
  }
  if (s_c < 7.348469228349534)
    // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,6) for g(y) := bₗ(-y²/2)/bₘₐₓ(-y²/2).
    // Branch III. Accuracy better than 7.49E-17 in perfect arithmetic.
    return (-9.3325115354837883291E-5 + s_c * (5.3118033972794648837E-4 + s_c * (7.4114855448345002595E-2 + s_c * (7.4039658186822817454E-2 + s_c * (3.9225177407687604785E-2 + s_c * (1.0022913378254090083E-2 + 1.7012579407246055469E-3 * s_c)))))) / (1 + s_c * (2.2217238132228132256 + s_c * (2.3441816707087403282 + s_c * (1.3912323646271141826 + s_c * (5.3231258443501838354E-1 + s_c * (1.1744005919716101572E-1 + 1.6195405895930935811E-2 * s_c))))));
#if defined( USE_RECIPROCAL_EVALUATION_IN_FAR_TAIL_OF_BL_OVER_BMAX )
  // x = -1/(2·y²), y = 1/√|2x|, output is f(x) = bₗ(x)/bₘₐₓ(x). Nonlinear-Remez optimized minimax rational function of order (6,4) for g(y) := ((f - Φ(-√(π/2)) )·4·exp(π/4)/y+√(π/2))/y .   f = Φ(-√(π/2)) + exp(-π//4)/4 · y · ( -√(π/2) + y · g ) .
  const double y = 1 / s_c, g = (4.6053948291412933458E-2 + y * (1.8614246236051740003 + y * (2.687739780481516225 + y * (2.1166111659498009072E1 + y * (7.7041195462942084146 + y * (1.364789754674820047E1 - 1.1914640506516351326E1 * y)))))) / (1 + y * (1.2314848469729795224 + y * (1.4972388415765707953E1 + y * (8.5172763098903471958 + 4.051186378326048589E1 * y))));
  // Branch IV. Accuracy better than 9.34E-17 in perfect arithmetic.
  // Φ(-√(π/2)) = 0.1050457027219686376856
  return 0.1050457027219686376856 + 0.1139845319414990591915 * y * (-1.253314137315500251208 + y * g);
#else
  // x = -1/(2·y²), y = 1/√|2x|, output is f(x) = bₗ(x)/bₘₐₓ(x). Remez-optimized minimax rational function of order (5,6) for g(y) := (f - Φ(-√(π/2))/y .   f = Φ(-√(π/2)) + y · g(y) .  Φ(-√(π/2)) = 0.1050457027219686376856.
  // The final solution is transformed back to a (6,6) rational function of s_c via full analytical simplification of Φ(-√(π/2)) + g(1/s_c) / s_c.
  // Branch IV. Accuracy better than 8.4E-17 in perfect arithmetic.
  return (1.4500072297240603183E-3 + s_c * (-1.5116692485011195757E-3 + s_c * (7.1682178310936334831E-2 + s_c * (3.921610857820463493E-2 + s_c * (2.9342405658628443931E-2 + s_c * (5.1832526171631521426E-3 + 1.6930208078421474854E-3 * s_c)))))) / (1 + s_c * (1.6176313502305414664 + s_c * (1.6823159175281531664 + s_c * (8.4878307567372222113E-1 + s_c * (3.7543742137375791321E-1 + s_c * (7.126137099644302999E-2 + 1.6116992546788676159E-2 * s_c))))));
#endif
#else
  if (s_c >= 2.6267851073127395) {
    if (s_c >= 7.348469228349534) {
#if defined( USE_RECIPROCAL_EVALUATION_IN_FAR_TAIL_OF_BL_OVER_BMAX )
      // x = -1/(2·y²), y = 1/√|2x|, output is f(x) = bₗ(x)/bₘₐₓ(x). Nonlinear-Remez optimized minimax rational function of order (6,4) for g(y) := ((f - Φ(-√(π/2)) )·4·exp(π/4)/y+√(π/2))/y .   f = Φ(-√(π/2)) + exp(-π//4)/4 · y · ( -√(π/2) + y · g ) .
      const double y = 1 / s_c, g = (4.6053948291412933458E-2 + y * (1.8614246236051740003 + y * (2.687739780481516225 + y * (2.1166111659498009072E1 + y * (7.7041195462942084146 + y * (1.364789754674820047E1 - 1.1914640506516351326E1 * y)))))) / (1 + y * (1.2314848469729795224 + y * (1.4972388415765707953E1 + y * (8.5172763098903471958 + 4.051186378326048589E1 * y))));
      // Branch IV. Accuracy better than 9.34E-17 in perfect arithmetic.
      // Φ(-√(π/2)) = 0.1050457027219686376856
      return 0.1050457027219686376856 + 0.1139845319414990591915 * y * (-1.253314137315500251208 + y * g);
#else
      // x = -1/(2·y²), y = 1/√|2x|, output is f(x) = bₗ(x)/bₘₐₓ(x). Remez-optimized minimax rational function of order (5,6) for g(y) := (f - Φ(-√(π/2))/y .   f = Φ(-√(π/2)) + y · g(y) .  Φ(-√(π/2)) = 0.1050457027219686376856.
      // The final solution is transformed back to a (6,6) rational function of s_c via full analytical simplification of Φ(-√(π/2)) + g(1/s_c) / s_c.
      // Branch IV. Accuracy better than 8.4E-17 in perfect arithmetic.
      return (1.4500072297240603183E-3 + s_c * (-1.5116692485011195757E-3 + s_c * (7.1682178310936334831E-2 + s_c * (3.921610857820463493E-2 + s_c * (2.9342405658628443931E-2 + s_c * (5.1832526171631521426E-3 + 1.6930208078421474854E-3 * s_c)))))) / (1 + s_c * (1.6176313502305414664 + s_c * (1.6823159175281531664 + s_c * (8.4878307567372222113E-1 + s_c * (3.7543742137375791321E-1 + s_c * (7.126137099644302999E-2 + 1.6116992546788676159E-2 * s_c))))));
#endif
    }
    // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,6) for g(y) := bₗ(-y²/2)/bₘₐₓ(-y²/2).
    // Branch III. Accuracy better than 7.49E-17 in perfect arithmetic.
    return (-9.3325115354837883291E-5 + s_c * (5.3118033972794648837E-4 + s_c * (7.4114855448345002595E-2 + s_c * (7.4039658186822817454E-2 + s_c * (3.9225177407687604785E-2 + s_c * (1.0022913378254090083E-2 + 1.7012579407246055469E-3 * s_c)))))) / (1 + s_c * (2.2217238132228132256 + s_c * (2.3441816707087403282 + s_c * (1.3912323646271141826 + s_c * (5.3231258443501838354E-1 + s_c * (1.1744005919716101572E-1 + 1.6195405895930935811E-2 * s_c))))));
  }
  if (s_c >= 0.7099295739719539)
    // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,6) for g(y) := bᵤ(-y²/2)/bₘₐₓ(-y²/2).
    // Branch II. Accuracy better than 8.77E-17 in perfect arithmetic.
    return (1.9795737927598581235E-9 + s_c * (-2.7081288564685588037E-8 + s_c * (7.5610142272549044609E-2 + s_c * (6.917130174466834016E-2 + s_c * (2.9537058950963019803E-2 + s_c * (6.5849252702302307774E-3 + 6.9711400639834715731E-4 * s_c)))))) / (1 + s_c * (2.1941448525586579756 + s_c * (2.1297103549995181357 + s_c * (1.1571483187179784072 + s_c * (3.7831622253060456794E-1 + s_c * (7.1714862448829349869E-2 + 6.6361975827861200167E-3 * s_c))))));
  // x = -y²/2, y = √|2x| = s_c, output is f(x) = bₗ(x)/bₘₐₓ(x).
  // For small |x|, i.e., small y, we have f ≈ (exp(-1/π)/4-Φ(-√(π/2))/2)·y² + exp(-1/π)/(3·√(2π))·y³ + O(y⁴).
  //   c₂ := (exp(-1/π)/4-Φ(-√(π/2))/2) =  0.07560996640296361767172
  //   c₃ := exp(-1/π)/(3·√(2π))        = -0.09672719281339436290858
  // Nonlinear-Remez optimized minimax rational function of order (5,4) for g(y) := ((f/y²-c₂)/y-c₃)/y .   f = y²·(c₂+y·(c₃+y·g)) .
  const double g = (8.0741072372882856924E-2 + s_c * (9.8078911786358897272E-2 + s_c * (3.9760631445677058375E-2 + s_c * (5.9716928459589189876E-3 + s_c * (-6.4036399341479799981E-6 + 4.5425102093616062245E-7 * s_c))))) / (1 + s_c * (1.8594977672287664353 + s_c * (1.3658801475711790419 + s_c * (4.6132707108655653215E-1 + 6.1254597049831720643E-2 * s_c))));
  // Branch I. Accuracy better than 7.43E-17 in perfect arithmetic.
  return (s_c * s_c) * (0.07560996640296361767172 + s_c * (s_c * g - 0.09672719281339436290858));
#endif
}

inline double bᵤ_over_bₘₐₓ(const double s_c /* = √|2x| */) {
  // Four branches.
  //                    │     I.          │        II.           │      III.        │        IV.     │
  //                    │    ←-→          │        ←-→           │      ←-→         │        ←-→     │
  //  In |x|:           [  0  ,  0.3                ,  1.6               ,  19                ,  ∞   )
  //  In s_c = √|2x| :  [  0  ,  0.7745966692414833 , 1.7888543819998317 , 6.164414002968976  ,  ∞   )
#if !defined(BU_BRANCHES_UNORDERED)
  if (s_c < 1.7888543819998317) {
    if (s_c < 0.7745966692414833) {
      // x = -y²/2, y = √|2x| = s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x).
      // For small |x|, i.e., small y, we have f ≈ 1-2·Φ(-√(π/2)) + (exp(-π/4)/4-Φ(-√(π/2))/2)·y² + O(y³).
      //   c₀ := 1-2·Φ(-√(π/2))             = 0.7899085945560627246288
      //   c₂ := (exp(-π/4)/4-Φ(-√(π/2))/2) = 0.0614616805805147403487
      // Nonlinear-Remez optimized minimax rational function of order (5,4) for g(y) := ((f - c₀)/y²-c₂)/y .   f = c₀ + y²·(c₂+y·g) .
      const double g = (-6.063099881233561706E-2 + s_c * (-8.1011946637120604985E-2 + s_c * (-4.2505564862438753828E-2 + s_c * (-8.9880000946868691788E-3 + s_c * (-7.5603072110443268356E-6 + 4.3879556621540147458E-7 * s_c))))) / (1 + s_c * (1.8400371530721828756 + s_c * (1.5709283443886143691 + s_c * (6.8913245453611400484E-1 + 1.4703173061720980923E-1 * s_c))));
      // Branch I. Accuracy better than 9.2E-17 in perfect arithmetic.
      return 0.7899085945560627246288 + (s_c * s_c) * (0.0614616805805147403487 + s_c * g);
    }
    // x = -y²/2, y = √|2x| = s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,5) for g(y) := bᵤ(-y²/2)/bₘₐₓ(-y²/2).
    // Branch II. Accuracy better than 8.4E-17 in perfect arithmetic.
    return (7.8990944435755287611E-1 + s_c * (-1.2655410534988972886 + s_c * (-2.8803040699221003256 + s_c * (-2.6936198689113258727 + s_c * (-1.1213067281643205754 + s_c * (-2.1277793801691629892E-1 + 5.1486445905299802703E-6 * s_c)))))) / (1 + s_c * (-1.6021222722060444448 + s_c * (-3.7242680976480704555 + s_c * (-3.2083117718907365085 + s_c * (-1.2922333835930958583 - 2.3762328334050001161E-1 * s_c)))));
  }
  if (s_c < 6.164414002968976) {
    // x = -y²/2, y = √|2x| = s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,6) for g(y) := bᵤ(-y²/2)/bₘₐₓ(-y²/2).
    // Branch III. Accuracy better than 7.7E-17 in perfect arithmetic.
    return (7.8990640048967596475E-1 + s_c * (1.5993699253596663678 + s_c * (1.6481729039140370242 + s_c * (9.8227188109869200166E-1 + s_c * (3.6313557966186936883E-1 + s_c * (7.8277036261179606301E-2 + 9.3404307364538726214E-3 * s_c)))))) / (1 + s_c * (2.0247407005640401446 + s_c * (2.0087454279103740489 + s_c * (1.1627561803056961973 + s_c * (4.2004672123723823581E-1 + s_c * (8.9130862793887234546E-2 + 1.0436767768858021717E-2 * s_c))))));
  }
#if defined( USE_RECIPROCAL_EVALUATION_IN_FAR_TAIL_OF_BU_OVER_BMAX )
  // x = -1/(2·y²), y = 1/√|2x| = 1/s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x). Nonlinear-Remez optimized minimax rational function of order (6,4) for g(y) := ((f - Φ(√(π/2)) )·4·exp(π/4)/y+√(π/2))/y + (1-π/4)² .   f = Φ(√(π/2)) + exp(-π//4)/4 · y · ( -√(π/2) + y · (g-(1-π/4)²) ) .
  const double y = 1 / s_c, g = (1.124338374937498095E-10 + y * (1.8047098623684767751 + y * (2.5589524137987524953 + y * (1.9974228630749728307E1 + y * (8.7923365534769864932 + y * (1.0442103559363093567E1 - 9.0656649638740373552E0 * y)))))) / (1 + y * (1.2936587104647790675 + y * (1.4215562205207917075E1 + y * (8.0353139924329651631 + 3.5848286453520043747E1 * y))));
  // Branch IV. Accuracy better than 8.9E-17 in perfect arithmetic.
  return 0.8949542972780313623144 + 0.1139845319414990591915 * y * (-1.253314137315500251208 + y * (g - 0.04605394827318829444583));
#else
  // x = -1/(2·y²), y = 1/√|2x| = 1/s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x). Remez-optimized minimax rational function of order (5,6) for g(y) := (f - Φ(√(π/2))/y .   f = Φ(√(π/2)) + y · g(y) .  Φ(√(π/2)) = 0.8949542972780313623144 .
  // The final solution is transformed back to a (6,6) rational function of s_c via full analytical simplification of Φ(√(π/2)) + g(1/s_c) / s_c.
  // Branch IV. Accuracy better than 3.9E-17 in perfect arithmetic.
  return (7.91133825948419359E-1 + s_c * (1.24653733210880042 + s_c * (1.32747426980537386 + s_c * (6.95009705717846778E-1 + s_c * (3.05965944268228457E-1 + s_c * (6.02200363391352887E-2 + 1.29050244454344842E-2 * s_c)))))) / (1 + s_c * (1.58117486714634672 + s_c * (1.60144713247629644 + s_c * (8.30040185836882436E-1 + s_c * (3.53071863813401531E-1 + s_c * (6.95901684131758475E-2 + 1.44197580643890011E-2 * s_c))))));
#endif
#else
  if (s_c > 1.7888543819998317) {
    if (s_c > 6.164414002968976) {
#if defined( USE_RECIPROCAL_EVALUATION_IN_FAR_TAIL_OF_BU_OVER_BMAX )
      // x = -1/(2·y²), y = 1/√|2x| = 1/s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x). Nonlinear-Remez optimized minimax rational function of order (6,4) for g(y) := ((f - Φ(√(π/2)) )·4·exp(π/4)/y+√(π/2))/y + (1-π/4)² .   f = Φ(√(π/2)) + exp(-π//4)/4 · y · ( -√(π/2) + y · (g-(1-π/4)²) ) .
      const double y = 1 / s_c, g = (1.124338374937498095E-10 + y * (1.8047098623684767751 + y * (2.5589524137987524953 + y * (1.9974228630749728307E1 + y * (8.7923365534769864932 + y * (1.0442103559363093567E1 - 9.0656649638740373552E0 * y)))))) / (1 + y * (1.2936587104647790675 + y * (1.4215562205207917075E1 + y * (8.0353139924329651631 + 3.5848286453520043747E1 * y))));
      // Branch IV. Accuracy better than 8.9E-17 in perfect arithmetic.
      return 0.8949542972780313623144 + 0.1139845319414990591915 * y * (-1.253314137315500251208 + y * (g - 0.04605394827318829444583));
#else
      // x = -1/(2·y²), y = 1/√|2x| = 1/s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x). Remez-optimized minimax rational function of order (5,6) for g(y) := (f - Φ(√(π/2))/y .   f = Φ(√(π/2)) + y · g(y) .  Φ(√(π/2)) = 0.8949542972780313623144 .
      // The final solution is transformed back to a (6,6) rational function of s_c via full analytical simplification of Φ(√(π/2)) + g(1/s_c) / s_c.
      // Branch IV. Accuracy better than 3.9E-17 in perfect arithmetic.
      return (7.91133825948419359E-1 + s_c * (1.24653733210880042 + s_c * (1.32747426980537386 + s_c * (6.95009705717846778E-1 + s_c * (3.05965944268228457E-1 + s_c * (6.02200363391352887E-2 + 1.29050244454344842E-2 * s_c)))))) / (1 + s_c * (1.58117486714634672 + s_c * (1.60144713247629644 + s_c * (8.30040185836882436E-1 + s_c * (3.53071863813401531E-1 + s_c * (6.95901684131758475E-2 + 1.44197580643890011E-2 * s_c))))));
#endif
    }
    // x = -y²/2, y = √|2x| = s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,6) for g(y) := bᵤ(-y²/2)/bₘₐₓ(-y²/2).
    // Branch III. Accuracy better than 7.7E-17 in perfect arithmetic.
    return (7.8990640048967596475E-1 + s_c * (1.5993699253596663678 + s_c * (1.6481729039140370242 + s_c * (9.8227188109869200166E-1 + s_c * (3.6313557966186936883E-1 + s_c * (7.8277036261179606301E-2 + 9.3404307364538726214E-3 * s_c)))))) / (1 + s_c * (2.0247407005640401446 + s_c * (2.0087454279103740489 + s_c * (1.1627561803056961973 + s_c * (4.2004672123723823581E-1 + s_c * (8.9130862793887234546E-2 + 1.0436767768858021717E-2 * s_c))))));
  }
  if (s_c >= 0.7745966692414833) {
    // x = -y²/2, y = √|2x| = s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x). Remez optimized minimax rational function of order (6,5) for g(y) := bᵤ(-y²/2)/bₘₐₓ(-y²/2).
    // Branch II. Accuracy better than 8.4E-17 in perfect arithmetic.
    return (7.8990944435755287611E-1 + s_c * (-1.2655410534988972886 + s_c * (-2.8803040699221003256 + s_c * (-2.6936198689113258727 + s_c * (-1.1213067281643205754 + s_c * (-2.1277793801691629892E-1 + 5.1486445905299802703E-6 * s_c)))))) / (1 + s_c * (-1.6021222722060444448 + s_c * (-3.7242680976480704555 + s_c * (-3.2083117718907365085 + s_c * (-1.2922333835930958583 - 2.3762328334050001161E-1 * s_c)))));
  }
  // x = -y²/2, y = √|2x| = s_c, output is f(x) = bᵤ(x)/bₘₐₓ(x).
  // For small |x|, i.e., small y, we have f ≈ 1-2·Φ(-√(π/2)) + (exp(-π/4)/4-Φ(-√(π/2))/2)·y² + O(y³).
  //   c₀ := 1-2·Φ(-√(π/2))             = 0.7899085945560627246288
  //   c₂ := (exp(-π/4)/4-Φ(-√(π/2))/2) = 0.0614616805805147403487
  // Remez-optimized minimax rational function of order (5,4) for g(y) := ((f - c₀)/y²-c₂)/y .   f = c₀ + y²·(c₂+y·g) .
  const double g = (-6.063099881233561706E-2 + s_c * (-8.1011946637120604985E-2 + s_c * (-4.2505564862438753828E-2 + s_c * (-8.9880000946868691788E-3 + s_c * (-7.5603072110443268356E-6 + 4.3879556621540147458E-7 * s_c))))) / (1 + s_c * (1.8400371530721828756 + s_c * (1.5709283443886143691 + s_c * (6.8913245453611400484E-1 + 1.4703173061720980923E-1 * s_c))));
  // Branch I. Accuracy better than 9.2E-17 in perfect arithmetic.
  return 0.7899085945560627246288 + (s_c * s_c) * (0.0614616805805147403487 + s_c * g);
#endif
}

#endif

#endif

#if defined( SPECIALISE_EXACT_ATM_CASE )

#if defined( USE_SPECIALISED_RATIONAL_FUNCTION_APPROXIMATION_FOR_ATM_IMPLIED_VOLATILITY )

#if !defined( USE_ALGORITHM_AS241 )
double inverse_norm_cdf_for_low_probabilities(double p);
#endif

// Specialisation for x = 0 where bₐₜₘ(s) = 1-2·Φ(-s/2) = 2·Φ(s/2)-1 = erf(s/√8).
double implied_normalised_volatility_atm(double 𝛽) {
  const double 𝛽ₘₐₓ = 0.6826894921370859; // = 1-erfc(1/√2)
  if (𝛽 <= 𝛽ₘₐₓ) {
    // Note that this branch domain corresponds exactly to innermost branch of Φ⁻¹() in algorithm "PJ-2024-Inverse-Normal" denoted as 'inverse_norm_cdfmhalf_for_midrange_probabilities()'
    // as also used in erfinv() in normaldistribution.cpp. The only difference is that the computation of s via erfinv() effectively does
    //   s = √8 · ( inverse_norm_cdfmhalf_for_midrange_probabilities(𝛽/2) / √2 ) .
    // For very small values of 𝛽, say DBL_MIN ≤ 𝛽 < 3·DBL_MIN or so, the innermost division by two in the term 𝛽/2 can already go below DBL_MIN
    // and thus becomes effectively 0 when we are compiling with -Ofast. For such small values, we should, however, have
    //   s = √(2π) · 𝛽
    //     ≈ 2.507 · 𝛽
    // i.e., s > 𝛽 by more than a factor of two. Expanding inverse_norm_cdfmhalf_for_midrange_probabilities(y) for very small y does indeed give
    //   inverse_norm_cdfmhalf_for_midrange_probabilities(y) ≈ √(2π) · y
    // which means in perfect arithmetic we'd indeed obtain
    //  √8 · erfinv(𝛽/2) = √8 · ( (√(2π)·𝛽/2) / √2 )
    //                   = √(2π) · 𝛽
    // Only that we obtain 0 for such small values of 𝛽 due to the intermediate underflow.
    // 
    // In practice, however, this only means that (for x=0) the minimum normalised Black price for which we can compute a positive implied
    // volatility and get the same normalised price back (i.e., we can exactly complete the full loop), is 6.67522e-308 instead of 4.45015e-308.
    // As for speed, I was unable to get a measurable difference between this specialised rational function approximation and the route via erfinv().
    // Given that there is no speed penalty, 6.67522e-308 vs. 4.45015e-308 as the minimum viable input price is an acceptable limitation to justify
    // the code simplification to go via the (otherwise reusable) erfinv() implementation.
    const double r = 𝛽ₘₐₓ * 𝛽ₘₐₓ - 𝛽 * 𝛽;
    // Remez - optimized minimax rational function approximation of order(6, 6) in r.
    // Accuracy better than 3.5E-17 in perfect arithmetic on β ∈ [0,1-erfc(1/√2)] which corresponds to s ∈ [0,2].
    return 𝛽 * ((2.92958954698308816 + r * (1.4014698674754995E1 + r * (2.44918990556468762E1 + r * (1.90763928424894996E1 + r * (6.43250149461895996 + r * (7.52328633671821543E-1 + 1.38781536163865582E-2 * r)))))) / (1 + r * (5.22443271807813073 + r * (1.02258209975070629E1 + r * (9.28187483709036392 + r * (3.9095549184069553 + r * (6.61214199809055912E-1 + 2.89411828874884851E-2 * r)))))));
  }
#if defined( USE_ALGORITHM_AS241 )
  return -2 * inverse_norm_cdf(0.5 * (1 - 𝛽));
#else
  // We fallback to the tail branch(es) of Φ⁻¹() in algorithm "PJ-2024-Inverse-Normal" [see file normaldistribution.cpp].
  return -2 * inverse_norm_cdf_for_low_probabilities(0.5 * (1 - 𝛽));
#endif
}
#else

// We can use the internal branches of Φ⁻¹(·) to implement erfinv() avoiding catastrophic subtractive cancellation for small arguments.
//double erfinv(double);

// Specialisation for x = 0 where bₐₜₘ(s) = 1-2·Φ(-s/2) = 2·Φ(s/2)-1 = erf(s/√8).
#define implied_normalised_volatility_atm(𝛽) ((2 * SQRT_TWO) * erfinv(𝛽))

#endif

#endif

// See http://en.wikipedia.org/wiki/Householder%27s_method for a detailed explanation of the third order Householder iteration.
//
// Given the objective function g(s) whose root x such that 0 = g(s) we seek, iterate
//
//     s_n+1  =  s_n  -  (g/g') · [ 1 - (g''/g')·(g/g') ] / [ 1 - (g/g')·( (g''/g') - (g'''/g')·(g/g')/6 ) ]
//
// Denoting  ν:=-(g/g'), h₂:=(g''/g'), and hh3:=(g'''/g'), this reads
//
//     s_n+1  =  s_n  +  ν · ( 1 + ν·h₂/2 ) / ( 1 + ν·( h₂ + ν·h₃/6 ) ).
//
//
// NOTE that this function requires that the input is for an out-of-the-money option, i.e. θ·x ≤ 0.
//
extern "C" double lets_be_rational(double 𝛽, double θx, int N) {
  assert(θx <= 0);
  if (𝛽 <= 0)
    return implied_volatility_output(0, 0 == 𝛽 ? 0 : VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC);
#if defined( POSITIVE_DENORMALISATION_CUTOFF )
  if (𝛽 < POSITIVE_DENORMALISATION_CUTOFF) // For positive but denormalised (a.k.a. 'subnormal') prices, we return 0 since it would be impossible to converge to full machine accuracy anyway.
    return implied_volatility_output(0, 0);
#endif
  const double bₘₐₓ = exp(0.5 * θx);
  if (𝛽 >= bₘₐₓ)
    return implied_volatility_output(0, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM);
  // Exactly at the forward is a not-so-uncommon real world use case.
#if defined( SPECIALISE_EXACT_ATM_CASE )
  //  bₐₜₘ(s) = 1-2·Φ(-s/2) = 2·Φ(s/2)-1 = erf(s/√8)
  if (0 == θx) return implied_normalised_volatility_atm(𝛽);
#endif
  int iterations = 0;
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT)
  int direction_reversal_count = 0;
  double ds_previous = 0;
#endif
  double f = -DBL_MAX, s = -DBL_MAX, ds = -DBL_MAX;
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT) || defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
  double s_left = DBL_MIN, s_right = DBL_MAX;
#endif
  // We have
  //  s_c  =  √(2|x|)   and  θ·x ≤ 0.
  //  b_c  =  b(θ·x,s_c)
  //       =  exp(θx/2)·Φ(θx/√(2|x|)+√(2|x|)/2) - exp(-θx/2)·Φ(θx/√(2|x|)-√(2|x|)/2)
  //       =  exp(-|x|/2)·Φ(-|x|/√(2|x|)+√(2|x|)/2) - exp(|x|/2)·Φ(-|x|/√(2|x|)-√(2|x|)/2)
  //       =  exp(-|x|/2)·Φ(0) - exp(|x|/2)·Φ(-√(2|x|))
  //       =  exp(-|x|/2) · [ 1/2 - exp(|x|)·Φ(-√(2|x|)) ]
  // ⟹ using Φ(z) = erfc(-z/√2)/2
  //  b_c  =  exp(-|x|/2) / 2 · [ 1 - exp(|x|)·erfc(√|x|) ]
  // ⟹ using erfc(y) = erfcx(y)·exp(-y²)
  //  b_c  =  exp(-|x|/2) / 2 · [ 1 - erfcx(√|x|) ]
  // Also:
  //  vega(x,s) = ∂b(x,s)/∂s = b'(s) = exp(-half·((x/s)²+(s/2)²) / √(2π)
  // So
  //  vega(x,s_c) = exp(-half·((x/√(2|x|))²+(√(2|x|)/2)²) / √(2π)
  //              = exp(-half·(|x|/2+|x|/2) / √(2π)
  //              = exp(-|x|/2) / √(2π)
  // Note: bₘₐₓ = exp(-|x|/2).
  const double sqrt_ax = sqrt(-θx), s_c = SQRT_TWO * sqrt_ax, ome = one_minus_erfcx(sqrt_ax), b_c = 0.5 * bₘₐₓ * ome;
  // Four branches.
  if (𝛽 < b_c) {   // LOWER HALF: s < s_c
    assert(θx < 0); // Even without ATM (x=0) specialisation, we cannot get here when x=0 since then b_c = 0. This assertion is merely a reminder and a safeguard against development errors.
    //  s_c  =  √(2|θx|)
    //  b_c  =  b(θx,s_c)  =  exp(θx/2)·Φ(-√|x|/√2+√|x|/√2) - exp(-θx/2)·Φ(-√|x|/√2-√|x|/√2)  =  half·exp(θx/2) - exp(-θx/2)·Φ(-√(2|x|))
    //  sₗ  :=   s_c - b(x,s_c)/b'(x,s_c)                         =   √(2|x|) - (exp(θx/2)·Φ(0) - exp(-θx/2)·Φ(-√(2|x|)))·exp((|x|/2+|x|/2)/2)·√(2π)
    //      =   √(2|x|) - [ half - exp(|x|)·Φ(-√(2|x|)) ] · √(2π)   =   √(2|x|) - √(π/2) · [ 1 - erfcx(√|x|) ]
    // For x ⟶ -∞, we get
    //   sₗ ⟶ √(2|x|) - √(π/2)
    // and for 
    //  bₗ  := b(θx,sₗ)
    // we get
    //   bₗ/bₘₐₓ ⟶ Φ(-√(π/2)) = 0.1050457027219686376856
    const double sₗ = s_c - SQRT_PI_OVER_TWO * ome; // = s_c - b_c / v_c;
    assert(sₗ > 0); // Even without ATM (x=0) specialisation, we are guaranteed to have a positive number here. This assertion is merely a reminder and a safeguard against development errors.
#if defined( USE_UNIVARIATE_RATIONAL_FUNCTION_APPROXIMATIONS_FOR_BL_AND_BU_OVER_BMAX )
    const double bₗ = bₗ_over_bₘₐₓ(s_c) * bₘₐₓ;
#else
    const double bₗ = normalised_black(θx, sₗ);
#endif
    if (𝛽 < bₗ) {   // LOWEST BRANCH:   s < sₗ
      {
        double f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2;
        compute_f_lower_map_and_first_two_derivatives(θx, sₗ, f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2);
        const double rₗₗ = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0., bₗ, 0., f_lower_map_l, 1., d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2, true);
        f = rational_cubic_interpolation(𝛽, 0., bₗ, 0., f_lower_map_l, 1., d_f_lower_map_l_d_beta, rₗₗ);
        if (!(f > 0)) { // This can happen due to roundoff truncation for extreme values such as |x|>500.
          // We switch to quadratic interpolation using f(0)≡0, f(bₗ), and f'(0)≡1 to specify the quadratic.
          const double t = 𝛽 / bₗ;
          f = (f_lower_map_l * t + bₗ * (1 - t)) * t;
        }
        s = inverse_f_lower_map(θx, f);
      }
      assert(s > 0); // Even without ATM (x=0) specialisation, we are guaranteed to have a positive number here. This assertion is merely a reminder and a safeguard against development errors.
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT) || defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
      s_right = sₗ;
#endif
      //
      // In this branch, which comprises the lowest segment, the objective function is
      //     g(s) = 1/ln(b(θx,s)) - 1/ln(𝛽)
      //          ≡ 1/ln(b(s)) - 1/ln(𝛽)
      // This makes
      //              g'                =   -b'/(b·ln(b)²)
      // using λ:=1/ln(b)
      //              g'                =   -b'/b·λ²
      //              ν      = -g/g'    =   (ln(𝛽)-ln(b))·ln(b)/ln(𝛽)·b/b'
      //                                =   (ln(𝛽)-1/λ)/(ln(𝛽)·λ) · b/b'     =   (λ-1/ln(𝛽))·b/(b'·λ²)
      //              h₂     = g''/g'   =   b''/b'  -  b'/b·(1+2/ln(b))
      //                                =   b''/b'  -  b'/b·(1+2·λ)
      //              h₃     = g'''/g'  =   b'''/b' +  2(b'/b)²·(1+3/ln(b)·(1+1/ln(b)))  -  3(b''/b)·(1+2/ln(b))
      //                                =   b'''/b' +  (b'/b)²·(2+6·λ·(1+λ))  -  (b''/b)·3·(1+2·λ)
      //              h₄     = g''''/g' =   b''''/b' - (b'/b)³·(6+λ·(22+λ·(36+λ·24))) + (b'/b)·(b''/b)·(12+36·λ·(1+λ)) - (b''/b)·(b''/b')·(3+6·λ)  - (b'''/b)·(4+8·λ)
      //                                =   b''''/b' - (b'/b)·[ (b'/b)²·(6+λ·(22+λ·(36+λ·24))) - (b''/b)·(12+36·λ·(1+λ)) ] - (b''/b)·(b''/b')·3·(1+2·λ)  - (b'''/b)·4·(1+2·λ)
      //                                =   b_h₄ - (b'/b) · [ (b'/b)²·(6+λ·(22+λ·(36+λ·24))) - (b''/b)·(12+36·λ·(1+λ)) ] - (b''/b)·b_h₂·3·(1+2·λ)  - (b'''/b)·4·(1+2·λ)
      // The Householder(3) iteration is
      //     s_n+1  =  s_n  +  ν · ( 1 + ν·h₂/2 ) / ( 1 + ν·( h₂ + ν·h₃/6 ) ).
      //
      const double ln_𝛽 = log(𝛽);
      for (; iterations<N && fabs(ds)>DBL_EPSILON * s; ++iterations) {
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT)
        if (ds * ds_previous < 0)
          ++direction_reversal_count;
        if (N > 3 && iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right))) {
          // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
          // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
#ifdef LOG_BINARY_NESTING
          if (direction_reversal_count > 2)
            printf("Intercepting excessive direction reversal in lowest branch.\n");
          else
            printf("Intercepting bracket boundary contact/violation in lowest branch.\n");
#endif
          s = 0.5 * (s_left + s_right);
          if (s_right - s_left <= DBL_EPSILON * s) break;
          direction_reversal_count = 0;
          ds = 0;
        }
        ds_previous = ds;
#endif
        assert(s > 0);
        // Structured bindings trigger a warning about the need for -std=c++17 with some older versions of g++. I don't like the tuple syntax: I very much prefer structured binding syntax, but I dislike warnings even more.
#if __cplusplus >= 201703L
        const auto [bx, ln_vega] = scaled_normalised_black_and_ln_vega(θx, s);
#else
        const auto bx_and_ln_vega = scaled_normalised_black_and_ln_vega(θx, s);
        const double bx = std::get<0>(bx_and_ln_vega), ln_vega = std::get<1>(bx_and_ln_vega);
#endif
        const double ln_b = log(bx) + ln_vega, bpob = 1 / bx;
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT)
        const double b = exp(ln_b), bp = bpob * b;
        if (b > 𝛽 && s < s_right) s_right = s; else if (b<𝛽 && s>s_left) s_left = s; // Tighten the bracket if applicable.
        if (!(b > 0 && bp > 0)) { // Numerical underflow. Switch to binary nesting for this iteration.
#ifdef LOG_BINARY_NESTING
          printf("Binary nesting in lowest branch: b=%g, b'=%g.\n", b, bp);
#endif
          ds = 0.5 * (s_left + s_right) - s;
        } else
#endif
        {
          const double h = θx / s, x²_over_s³ = h * h / s, b_h₂ = x²_over_s³ - s / 4, ν = (ln_𝛽 - ln_b) * ln_b / ln_𝛽 / bpob, λ = 1 / ln_b, otλ = 1 + 2 * λ, h₂ = b_h₂ - bpob * otλ, c = 3 * (x²_over_s³ / s) /* = (h/s)² */;
          const double b_h₃ = b_h₂ * b_h₂ - c - 0.25, sq_bpob = bpob * bpob, bppob = b_h₂ * bpob, μ = 6 * λ * (1 + λ), h₃ = b_h₃ + sq_bpob * (2 + μ) - bppob * 3 * otλ;
          //
          // Introduced on 2023-12-14: for very large moneyness ratios [of no commercial relevance: exp(-190) = 3.05E-83], with exactly two Householder(3) iterations, I noticed that there is systematically a
          // residual inaccuracy in this branch [0 < b < bₗ] higher than the theoretically attainable one given by (|b(s)/(s·b'(s))|+1)·ε where ε is DBL_EPSILON and b(s) is the normalised Black function.
          // This residual inaccuracy disappears when we use two Householder(4) [5th order accuracy] iterations instead. Tests show that the initial guess is always close enough for the method to be contractive.
          // See further down in the description of the BlackAccuracyFactor() for a derivation of the formula (|b(s)/(s·b'(s))|+1)·ε. In this branch, we find that (|b(s)/(s·b'(s))|+1) is numerically equal to 1,
          // and thus the theoretically attainable relative accuracy is DBL_EPSILON.
          // 
          if (θx < -190) {
            //  b_h₄   =    b''''/b'   =   b_h₂·(b_h₃-half) - (b_h₂-2/s)·6·x²/s⁴
            //   h₄    =    b_h₄ - (b'/b) · [ (b'/b)²·(6+λ·(22+λ·(36+λ·24))) - (b''/b)·(12+36·λ·(1+λ)) ] - (b''/b)·b_h₂·3·(1+2·λ)  - (b'''/b)·4·(1+2·λ)    with    λ=1/ln(b)
            ds = ν * householder_factor(ν, h₂, h₃, (b_h₂ * (b_h₃ - 0.5) - (b_h₂ - 2 / s) * 2 * c) - bpob * (sq_bpob * (6 + λ * (22 + λ * (36 + λ * 24))) - bppob * (12 + 6 * μ)) - bppob * b_h₂ * 3 * otλ - b_h₃ * bpob * 4 * otλ);
          } else
            ds = ν * householder_factor(ν, h₂, h₃);
        }
#if defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
        // Never leave the branch (or bracket)
        s = take_step(s_left, s_right, s, ds);
#else
        s += ds;
#endif
        assert(s > 0);
      }
      return implied_volatility_output(iterations, s);
    } else {  // Lower middle: sₗ ≤ s < s_c
      const double inv_v_c = SQRT_TWO_PI / bₘₐₓ;  // v_c = bₘₐₓ * (1 / SQRT_TWO_PI)
      const double inv_vₗ = inv_normalised_vega(θx, sₗ), rₗₘ = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(bₗ, b_c, sₗ, s_c, inv_vₗ, inv_v_c, 0.0, false);
      s = rational_cubic_interpolation(𝛽, bₗ, b_c, sₗ, s_c, inv_vₗ, inv_v_c, rₗₘ);
      assert(s > 0); // Even without ATM (x=0) specialisation, we are guaranteed to have a positive number here. This assertion is merely a reminder and a safeguard against development errors.
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT) || defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
      s_left = sₗ;
      s_right = s_c;
#endif
    }
  } else {         // UPPER HALF:  s_c ≤ s
    //  sᵤ  :=   s_c + (bₘₐₓ-b(x,s_c))/b'(x,s_c)
    //       =   √(2|x|) + (exp(θx/2) - exp(θx/2)·Φ(0) + exp(-θx/2)·Φ(-√(2|x|))) · exp((|θx|/2+|θx|/2)/2)·√(2π)
    //       =   √(2|x|) + (exp(θx/2)/2 + exp(-θx/2)·Φ(-√(2|x|))) · exp(-θx/2)·√(2π)
    //       =   √(2|x|) + [ half + exp(-θx)·Φ(-√(2|x|)) ] · √(2π)
    //       =   √(2|x|) + √(π/2) · [ 1 + erfcx(√|x|) ]
    //       =   sₗ + √(2π)
    // For x ⟶ -∞, we get
    //   sᵤ ⟶ √(2|x|) + √(π/2)
    // and for
    //  bᵤ  := b(θx,sᵤ)
    // we get
    //   bᵤ/bₘₐₓ ⟶ Φ(√(π/2)) = 0.8949542972780313623144
    //
#if defined( DO_NOT_SIMPLIFY_SU )
    const double inv_v_c = SQRT_TWO_PI / bₘₐₓ;  // v_c = bₘₐₓ * (1 / SQRT_TWO_PI)
    const double sᵤ = v_c > DBL_MIN ? s_c + (bₘₐₓ - b_c) * inv_v_c : s_c;
#else
    const double sᵤ = s_c + SQRT_PI_OVER_TWO * (2 - ome); // = s_c + (bₘₐₓ - b_c) / v_c    ---     ome  =  1 - erfcx(√|x|)
#endif
    assert(sᵤ > 0);
#if defined( USE_UNIVARIATE_RATIONAL_FUNCTION_APPROXIMATIONS_FOR_BL_AND_BU_OVER_BMAX )
    const double bᵤ = bᵤ_over_bₘₐₓ(s_c) * bₘₐₓ;
#else
    const double bᵤ = normalised_black(θx, sᵤ);
#endif
    if (𝛽 <= bᵤ) { // UPPER MIDDLE:  s_c ≤ s ≤ sᵤ
#if !defined( DO_NOT_SIMPLIFY_SU )
      const double inv_v_c = SQRT_TWO_PI / bₘₐₓ;  // v_c = bₘₐₓ * (1 / SQRT_TWO_PI)
#endif
      const double inv_vᵤ = inv_normalised_vega(θx, sᵤ), rᵤₘ = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_c, bᵤ, s_c, sᵤ, inv_v_c, inv_vᵤ, 0.0, false);
      s = rational_cubic_interpolation(𝛽, b_c, bᵤ, s_c, sᵤ, inv_v_c, inv_vᵤ, rᵤₘ);
      assert(s > 0);
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT) || defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
      s_left = s_c;
      s_right = sᵤ;
#endif
    } else {       // HIGHEST BRANCH:  sᵤ < s  and  𝛽 > bₘₐₓ/2
      // The target value 𝛽 ϵ [bᵤ,bₘₐₓ).
      double f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2;
      compute_f_upper_map_and_first_two_derivatives(θx, sᵤ, f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2);
      if (d2_f_upper_map_h_d_beta2 > -SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2 < SQRT_DBL_MAX) {
        const double rᵤᵤ = convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(bᵤ, bₘₐₓ, f_upper_map_h, 0., d_f_upper_map_h_d_beta, -0.5, d2_f_upper_map_h_d_beta2, true);
        f = rational_cubic_interpolation(𝛽, bᵤ, bₘₐₓ, f_upper_map_h, 0., d_f_upper_map_h_d_beta, -0.5, rᵤᵤ);
      }
      if (f <= 0) {
        const double h = bₘₐₓ - bᵤ, t = (𝛽 - bᵤ) / h;
        f = (f_upper_map_h * (1 - t) + 0.5 * h * t) * (1 - t); // We switch to quadratic interpolation using f(bᵤ), f(bₘₐₓ)≡0, and f'(bₘₐₓ)≡-1/2 to specify the quadratic.
      }
      s = inverse_f_upper_map(f);
      assert(s > 0);
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT) || defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
      s_left = sᵤ;
#endif
      if (𝛽 > 0.5 * bₘₐₓ) { // Else we better drop through and let the objective function be g(s) = b(θx,s)-𝛽. 
        //
        // In this branch, which comprises the upper segment, the objective function is
        //     g(s) = ln(bₘₐₓ-𝛽)-ln(bₘₐₓ-b(θx,s))
        //          ≡ ln((bₘₐₓ-𝛽)/(bₘₐₓ-b(s)))
        // This makes
        //              g'         =   b'/(bₘₐₓ-b)
        // 
        // from here on (see further down), using
        // 
        //           b̄(θx,s)      :=   bₘₐₓ   -  b(θx,s)
        // and
        //           β̄            :=   bₘₐₓ   -  𝛽
        //
        // we get for ν=-g/g', h₂=g''/g', h₃=g'''/g', h₄=g''''/g' :
        //
        //         ν   =  -g/g'     =   ln(b̄/β̄)·b̄/b'
        //         h₂  =  g''/g'    =   b''/b'   +  b'/b̄
        //                          =   b''/b'   +  g'
        //         h₃  =  g'''/g'   =   b'''/b'  +  g'·(2g'+3b''/b')
        //         h₄  =  g''''/g'  =   b''''/b' +  g'²·6·(g'+2b''/b') + 3·(b''/b')²·g' + 4·(b'''/b')·g'
        //             =  g''''/g'  =   b''''/b' +  g' · ( 6·g'·(g'+2b''/b') + 3·(b''/b')² + 4·(b'''/b') )
        // 
        // and the iteration is
        //     s_n+1  =  s_n  +  ν · ( 1 + ν·h₂/2 ) / ( 1 + ν·( h₂ + ν·h₃/6 ) ).
        //
        const double 𝛽_bar = bₘₐₓ - 𝛽;
        for (; iterations<N && fabs(ds)>DBL_EPSILON * s; ++iterations) {
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT)
          if (ds * ds_previous < 0)
            ++direction_reversal_count;
          if (N > 3 && iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right))) {
            // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
            // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
#ifdef LOG_BINARY_NESTING
            if (direction_reversal_count > 2)
              printf("Intercepting excessive direction reversal in highest branch.\n");
            else
              printf("Intercepting bracket boundary contact/violation in highest branch.\n");
#endif
            s = 0.5 * (s_left + s_right);
            if (s_right - s_left <= DBL_EPSILON * s) break;
            direction_reversal_count = 0;
            ds = 0;
          }
          ds_previous = ds;
          // See below as to the reason behind and the derivation of the formula for b̄(θx,s).
          const double h = θx / s, t = s / 2, gp /* = bp / b_bar */ = (2 / SQRT_TWO_PI) / (erfcx_cody((t + h) * (1 / SQRT_TWO)) + erfcx_cody((t - h) * (1 / SQRT_TWO))), bp = normalised_vega(θx, s), b_bar = bp / gp;
          // b > 𝛽  <=>  b̄ < β̄ and vice versa.
          if (b_bar < 𝛽_bar && s < s_right) s_right = s; else if (b_bar > 𝛽_bar && s > s_left) s_left = s; // Tighten the bracket if applicable.
          if (!(b_bar > DBL_MIN && bp > DBL_MIN)) { // Numerical over-/underflow. Switch to binary nesting for this iteration.
            // NOTE (2023-12-12): since the switch to the direct computation of b̄(θx,s) without any subtractive cancellation, I have no longer seen this branch entered into.
#ifdef LOG_BINARY_NESTING
            printf("Binary nesting in highest branch.\n");
#endif
            ds = 0.5 * (s_left + s_right) - s;
          } else
#else
          // See below as to the reason behind and the derivation of the formula for b̄(θx,s).
          const double h = θx / s, t = s / 2, gp /* = bp / b_bar */ = (2 / SQRT_TWO_PI) / (erfcx_cody((t + h) * (1 / SQRT_TWO)) + erfcx_cody((t - h) * (1 / SQRT_TWO))), b_bar = normalised_vega(θx, s) / gp;
#endif
          {
            //
            // Introduced on 2023-12-10
            //
            //    b̄(θx,s)           :=   bₘₐₓ   -  b(θx,s)
            //    b̄(θx,s)            =   exp(θx/2)  -  [ exp(θx/2)·Φ(θx/s+s/2) - exp(-θx/2)·Φ(θx/s-s/2) ]                     |     using  bₘₐₓ = exp(θx/2)
            //    b̄(θx,s)            =   exp(θx/2)·Φ(-θx/s-s/2)  +  exp(-θx/2)·Φ(θx/s-s/2)                       (1)     |     using  1-Φ(z) = Φ(-z)
            // 
            // Note: no subtractive cancellation!
            // 
            //    b̄(θx,s)            =   half · [ erfc((s/2+θx/s)/√2)·bₘₐₓ + erfc((s/2-θx/s)/√2)/bₘₐₓ ]   (2)     |     using  Φ(z) = erfc(-z/√2)/2
            // 
            // In this upper segment, b > bᵤ = b(su), sᵤ = sc + b̄(sc)/b'(sc)   with sc = √(2|x|) and dropping the dependency on x, we benefit from
            // the [re-]evaluation of (bₘₐₓ-b) via formula (1) or (2) above.
            // 
            // ················································································································································································································
            // To see this, consider that bᵤ(-|x|)/bₘₐₓ is monotonically increasing on |x| ∈ [0,∞). The bounds can be readily computed as follows, assuming w.l.o.g. that θx<0:
            //
            //      sᵤ(θx)   =   sc + b̄(θx,sc)/b'(sc)
            //               =   sc + [ exp(θx/2)·Φ(-θx/sc-sc/2)  +  exp(-θx/2)·Φ(θx/sc-sc/2) ] · exp(half·((θx/sc)²+(sc/2)²) · √(2π)               |   θx = - √(|x|·|x|), θx/sc = - √(|x|/2), sc/2 = +√(|x|/2), x/sc+sc/2 = 0, θx/sc-sc/2 = -sc
            //               =   √(2|x|) + [ half·exp(θx/2)  +  exp(-θx/2)·Φ(-√(2|x|)) ] · exp(half·|x|) · √(2π)
            //               =   √(2|x|)  +  √(π/2)  +  e⁻ˣ·Φ(-√(2|x|)) · √(2π)
            //
            // Limiting cases:
            //
            //      sᵤ(0)    =   √(2π)
            //      bᵤ(0)    =   b(0,sᵤ(0))   =   Φ(sᵤ(0)/2)  -  Φ(-sᵤ(0)/2)
            //               =   Φ(√(π/2))  -  Φ(-√(π/2))
            //               =   0.7899085945560624
            // 
            //  θx ⟶ -∞:
            // 
            //      sᵤ(θx)   =   √(2|x|)  +  √(π/2)  +  e⁻ˣ·Φ(-√(2|x|)) · √(2π)
            //      sᵤ(θx)   ≈   √(2|x|)  +  √(π/2)  +  e⁻ˣ·φ(-√(2|x|)) · √(2π) / √(2|x|)                                    | Abramowitz & Stegun 26.2.12
            //               =   √(2|x|)  +  √(π/2)  +  e⁻ˣ·exp(-|x|) / √(2|x|)
            //               =   √(2|x|)  +  √(π/2)  +  1 / √(2|x|)
            //               =   √(2|x|) · ( 1  +  √(π/(4|x|)) +  1/(4|x|) )
            //    θx/sᵤ(θx)  =   -√(|x|/2) · ( 1  -  √(π/(4|x|)) +  O(1/|x|) )
            //      bᵤ(θx)   =   b(θx,sᵤ(θx))
            //               ≈   exp(θx/2)·Φ( -√(|x|/2)·(1-√(π/(4|x|))) + √(|x|/2)·(1+√(π/(4|x|))) )    exp(-θx/2)·Φ( -√(|x|/2)·(1-√(π/(4|x|))) - √(|x|/2)·(1+√(π/(4|x|))) )
            //               =   exp(θx/2)·Φ( √(π/2) )  -  exp(-θx/2)·Φ( -√(2|x|) )
            //               =   exp(θx/2)·Φ( √(π/2) )  -  1/√(4π|x|)                                                             | Abramowitz & Stegun 26.2.12
            // 
            // With bₘₐₓ(θx) = exp(θx/2), this means
            // 
            //    bᵤ(0)/bₘₐₓ(0)   =   Φ(√(π/2))  -  Φ(-√(π/2))   =   0.7899085945560624
            // 
            // and
            //      lim  bᵤ(x)/bₘₐₓ(x)      =      Φ(√(π/2))     =   0.8949542972780312
            //       x→-∞
            // 
            // In other words, on s ∈ [sᵤ,∞), where b  ∈ [bᵤ,bₘₐₓ)  [ 0 < bᵤ ≤ b < bₘₐₓ ], we always have b ≥ bᵤ > ¾·bₘₐₓ and thus b̄ = bₘₐₓ-b < bₘₐₓ/4.
            // And whenever any f̄ = fₘₐₓ-f is less than f/4 we incur less roundoff error in f̄ if we can compute f̄ directly without subtractive cancellation. □ (q.e.d.)
            // ················································································································································································································
            //
            // Continuing from equation (2), using erfc(z) = erfcx(z)·exp(-z²), we get
            // 
            //    b̄(θx,s)   =   half · ( erfcx((s/2+θx/s)/√2)·exp(-half((θx/s)²+θx+(s/2)²))·exp(θx/2) + erfcx((s/2-θx/s)/√2)·exp(-half((θx/s)²-θx+(s/2)²))·exp(-θx/2) )
            //
            //              =   half · ( erfcx((s/2+θx/s)/√2) + erfcx((s/2-θx/s)/√2) ) · exp(-half((θx/s)²+(s/2)²))
            //
            //              =   √(π/2) · [ erfcx((s/2+θx/s)/√2) + erfcx((s/2-θx/s)/√2) ] · ∂b(θx,s)/∂s
            //
            // Ergo, ∂b(θx,s)/∂s / b̄(θx,s)  =  √(2/π) / ( erfcx((s/2+θx/s)/√2) + erfcx((s/2-θx/s)/√2) )
            //
            const double g = log(𝛽_bar / b_bar), x²_over_s³ = (h * h) / s, b_h₂ = x²_over_s³ - s / 4, c = 3 * (x²_over_s³ / s) /* = (h/s)² */, b_h₃ = b_h₂ * b_h₂ - c - 0.25;
            const double ν = -g / gp, h₂ = b_h₂ + gp, h₃ = b_h₃ + gp * (2 * gp + 3 * b_h₂);
            //
            // Introduced on 2023-12-14: for very large moneyness ratios [of no commercial relevance: exp(-580) = 1.286E-252], with exactly two Householder(3) iterations, I noticed that there is systematically a
            // residual inaccuracy in this branch (bᵤ ≤ b < bₘₐₓ = exp(θx/2)) higher than the theoretically attainable one given by (|b(s)/(s·b'(s))|+1)·ε where ε is DBL_EPSILON and b(s) is the normalised Black function.
            // This residual inaccuracy disappears when we use two Householder(4) [5th order accuracy] iterations instead. Tests show that the initial guess is always close enough for the method to be contractive.
            // See further down in the description of the BlackAccuracyFactor() for a derivation of the formula (|b(s)/(s·b'(s))|+1)·ε.
            // 
            if (θx < -580) {
              //  b_h₄   =    b''''/b'   =   b_h₂·(b_h₃-half) - (b_h₂-2/s)·6·x²/s⁴
              //   h₄    =    b''''/b' +  g' · ( 6·g'·(g'+2b''/b') + 3·(b''/b')² + 4·(b'''/b') )   =   b_h₄ +  g' · ( 6·g'·(g'+2·b_h₂) + 3·b_h₂² + 4·b_h₃ )
              ds = ν * householder_factor(ν, h₂, h₃, (b_h₂ * (b_h₃ - 0.5) - (b_h₂ - 2 / s) * 2 * c) + gp * (6 * gp * (gp + 2 * b_h₂) + 3 * b_h₂ * b_h₂ + 4 * b_h₃));
            } else
              ds = ν * householder_factor(ν, h₂, h₃);
          }
#if defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
          // Never leave the branch (or bracket)
          s = take_step(s_left, s_right, s, ds);
#else
          s += ds;
#endif
          assert(s > 0);
        }
        return implied_volatility_output(iterations, s);
      }
    }
  }
  //
  // MIDDLE BRANCHES (ITERATION):  sₗ ≤ s  and ( s < sᵤ or  𝛽 ≤ bₘₐₓ/2 )
  // 
  // In this branch, which comprises the two middle segments, the objective function is g(s) = b(θx,s)-𝛽, or g(s) = b(s) - 𝛽, for short.
  // This makes
  //              ν    =   -g/g'     =  -(b-𝛽)/b'
  //              h₂   =   g''/g'    =    b''/b'      =   x²/s³-s/4
  //              h₃   =   g'''/g'   =    b'''/b'     =   h₂² - 3·x²/s⁴ - 1/4
  //              h₄   =   g''''/g'  =    b''''/b'    =   h₂·(h₃-half)-(h₂-2/s)·6·x²/s⁴     [ not actually used in this branch ]
  // 
  // and the iteration is
  //     s_n+1  =  s_n  +  ν · ( 1 + ν·h₂/2 ) / ( 1 + ν·( h₂ + ν·h₃/6 ) ).
  //
  for (; iterations<N && fabs(ds)>DBL_EPSILON * s; ++iterations) {
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT)
    if (ds * ds_previous < 0)
      ++direction_reversal_count;
    if (N > 3 && iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right))) {
      // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
      // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
#ifdef LOG_BINARY_NESTING
      if (direction_reversal_count > 2)
        printf("Intercepting excessive direction reversal in highest branch.\n");
      else
        printf("Intercepting bracket boundary contact/violation in highest branch.\n");
#endif
      s = 0.5 * (s_left + s_right);
      if (s_right - s_left <= DBL_EPSILON * s) break;
      direction_reversal_count = 0;
      ds = 0;
    }
    ds_previous = ds;
#endif
    assert(s > 0);
    const double b = normalised_black(θx, s), inv_bp = inv_normalised_vega(θx, s), ν = (𝛽 - b) * inv_bp, h = θx / s, x²_over_s³ = (h * h) / s, h₂ = x²_over_s³ - s * 0.25, h₃ = h₂ * h₂ - 3 * (x²_over_s³ / s) /* = (h/s)² */ - 0.25;
#if defined (ENABLE_CHANGING_THE_MAXIMUM_ITERATION_COUNT)
    if (b > 𝛽 && s < s_right) s_right = s; else if (b<𝛽 && s>s_left) s_left = s; // Tighten the bracket if applicable.
#endif
    ds = ν * householder_factor(ν, h₂, h₃);
#if defined( SAFEGUARD_STEP_AGAINST_BRACKETS )
    // Never leave the branch (or bracket)
    s = take_step(s_left, s_right, s, ds);
#else
    s += ds;
#endif
    assert(s > 0);
  }
  return implied_volatility_output(iterations, s);
}

double NormalisedBlack(double x, double s, double θ /* θ=±1 */) {
  // Specialisation for x = 0 where b(s) = 1-2·Φ(-s/2) = erf(s/√8).
  if (0 == x) return erf_cody((0.5 / SQRT_TWO) * s);
  return normalised_intrinsic(θ < 0 ? -x : x) + (s <= 0 ? 0 : normalised_black(-fabs(x), s)); /* Reciprocal-strike call-put equivalence */
}

double Black(double F, double K, double sigma, double T, double θ /* θ=±1 */) {
  const double s = sigma * sqrt(T);
  // Specialisation for x = 0 where b(s) = 1-2·Φ(-s/2) = erf(s/√8).
  if (K == F) return F * erf_cody((0.5 / SQRT_TWO) * s);
  // Map in-the-money to out-of-the-money
  return std::max(θ < 0 ? K - F : F - K, 0.0) + (s <= 0 ? 0 : (sqrt(F) * sqrt(K)) * normalised_black(-fabs(log(F / K)), s));
}

//    b̄(x,s,θ)          :=   bₘₐₓ(x,θ)   -  b(x,s,θ)
//                       =   exp(θ·x/2)  -  θ·[ exp(θx/2)·Φ(θ·(x/s+s/2)) - exp(-θx/2)·Φ(θ·(x/s-s/2)) ]                |     using  bₘₐₓ = exp(θ·x/2)
//
//                           ⎧   exp(θx/2)·[1-Φ(x/s+s/2)] + exp(-θx/2)·Φ(x/s-s/2)                                     |     when θ = +1
//                       =   ⎨ 
//                           ⎩   exp(θx/2)·Φ(-x/s-s/2)    + exp(-θx/2)·[1-Φ(-x/s+s/2)]                                |     when θ = -1
// 
//                       =   exp(θx/2)·Φ(-x/s-s/2) + exp(-θx/2)·Φ(x/s-s/2)                                            |     for both θ = ±1
// 
double ComplementaryNormalisedBlack(double x, double s) { return complementary_normalised_black(x / s, s / 2); }

double ImpliedBlackVolatility(double price, double F, double K, double T, double θ /* θ=±1 */) {
  if (price >= (θ < 0 ? K : F))
    return implied_volatility_output(0, VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM);
  const double μ = θ < 0 ? K - F : F - K; // Map in-the-money to out-of-the-money
  return lets_be_rational((μ > 0 ? price - μ : price) / (sqrt(F) * sqrt(K)), -fabs(log(F / K)), IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS) / sqrt(T);
}

double NormalisedImpliedBlackVolatility(double 𝛽, double x, double θ /* q=±1 */) {
  // Map in-the-money to out-of-the-money
  return lets_be_rational(𝛽 - normalised_intrinsic(θ < 0 ? -x : x), -fabs(x), IMPLIED_VOLATILITY_MAXIMUM_ITERATIONS);
}

double NormalisedVega(double x, double s) {
  const double ax = fabs(x);
  if (ax <= 0)
    return (1 / SQRT_TWO_PI) * exp(-0.125 * s * s);
  if (s <= 0 || s <= ax * SQRT_DBL_MIN)
    return 0;
  return normalised_vega(x, s);
}

double Vega(double F, double K, double sigma, double T) { return (sqrt(F) * sqrt(K)) * NormalisedVega(log(F / K), sigma * sqrt(T)) * sqrt(T); }

double Volga(double F, double K, double sigma, double T) { return (sqrt(F) * sqrt(K)) * normalised_volga(log(F / K), sigma * sqrt(T)) * T; }

double NormalisedVolga(double x, double s) { return normalised_volga(x, s); }

double DblEpsilon() { return DBL_EPSILON; }

double DblMin() { return DBL_MIN; }

double DblMax() { return DBL_MAX; }

// Floating point numbers have finite precision. The resolution limit is defined as the smallest positive number ε such that 1 and 1+ε still have distinct representations
// in the respectively used floating point model. For standard IEEE 754 double precision (64 bit, 53 bit mantissa), that's about 2.22E-16 and defined as DBL_EPSILON in C/C++.
// We thus have to always assume that any input x into any function f() comes as a number that is representative for some number in the range (x-ε·x,x+ε·x). We will denote
// the concept of 'some number in the range (x-δx,x+δx)' as x±δx.
// 
// Error propagation in generic function evaluation.
// =================================================
//    Given an input number x with associated absolute precision δx, the evaluation of a function f(x) will incur both the uncertainty in the input as well as the finite
// precision of the representation of the result of the evaluation. In other words, instead of y = f(x), by propagation as well as incurred evaluation imprecision, we have
//    y±δy = f(x±δx)·(1±ε)
// which is to lowest order
//         = f(x) ± f'(x)·δx ± f(x)·ε
// Given that the two uncertainty terms on the right hand side can accumulate, net, using y = f(x) as the target (infinite precision) result, this means 
//    |δy| = |f'(x)·δx| + |f(x)·ε|
// which brings us to the *relative* precision of y as
//    |δy/y| = |f'(x)·δx/f(x)| + ε
// IF the input precision on x was |δx| = |x|·ε, we arrive at
//    |δy/y| = ( |x·f'(x)/f(x)| + 1 ) · ε.
// 
// Error propagation in inverse function evaluation.
// =================================================
// Given a function g(y) that is the inverse of another function f(·) such that g(f(x)) = x [in perfect arithmetic], we obtain from y = f(x) and x = g(y) via the same
// logic as above
//    |δx/x|  =  ( |y·g'(y)/g(y)| + 1 ) · ε   =  ( |f(x)/(x·f'(x))| + 1 ) · ε.
// In other words, if the evaluation of y := f(x) incurs a precision deterioration given by the multiplicative factor ( 1 + γ ) [where γ>0], i.e., if the input (relative) precision in x
// was ε and the output (relative) precision in y was (1+γ)·ε, then the evaluation of x := f⁻¹(y) results in the (relative) precision worsening from η := |δy/y| to |δx/x| = (1+1/γ)·η.
// 
// Error propagation in inverse functions with offset.
// ===================================================
// Consider a function f(x) limited from above by its x→∞ asymptotic value fₘₐₓ. Consider that we have, at least in theory, access to the [infinite precision] complementary
// function f̄(x) := fₘₐₓ - f(x). Naturally, we have f̄(x)→0 for x → ∞ , though an inverse of f̄(x) can be very well defined numerically just as the inverse of exp(-x) is well
// defined as x = -ln(y) for y = exp(-x). Having clarified the above, we now in fact seek the inverse of f(x) which obviously satisfies f(x) = fₘₐₓ - f̄(x). In support of
// finding the inverse of y = f(x), we define ȳ := fₘₐₓ - y and the inverse of f̄ as ̄ḡ(·) = f̄⁻¹(·) which is evaluated from any input value y as x := ḡ(fₘₐₓ-y). Note that we
// need to pay attention that the evaluation ȳ = ȳ(y) = fₘₐₓ - y incurs itself from any input value y±δy the error propagation
//     ȳ±δȳ  =  (fₘₐₓ-(y±δy))·(1±ε)  =  (fₘₐₓ-y) ± δy ± (fₘₐₓ-y)·ε ± δy·ε.
// Setting |δy| = |y|·ε since y is an input value, and using ȳ=fₘₐₓ-y gives us (recall that we must alas always accumulate errors in absolute value) to lowest order in ε
//     |δȳ|  =  (|y|+|ȳ|)·ε ,
// i.e.,
//     |δȳ/ȳ|  =  (1+|y/ȳ|)·ε .
// We now compute the error propagation of the inverse function evaluation ḡ(ȳ) with ȳ=(fₘₐₓ-y) from the input y±δy with |δy| = |y|·ε to lowest order in ε as follows:
//     x±δx  =  ḡ(ȳ±δȳ)·(1±ε)
//           =  ḡ(ȳ±(|y|+|ȳ|)·ε)·(1±ε)
//           =  ḡ(ȳ) ± ḡ'(ȳ)·(|y|+|ȳ|)·ε ± ḡ(ȳ)·ε
// Since ḡ(ȳ)=x, this yields
//     |δx|  =  [ |ḡ'(ȳ)|·(|y|+|ȳ|) + |x| ] · ε.
// Using ḡ'(ȳ) = 1/f̄'(x) and |f̄'(x)| = |f'(x)|, we continue
//     |δx|  =  [ (|y|+|ȳ|)/|f'(x)| + |x| ] · ε
// whence
//   |δx/x|  =  [ |f(x)|/|x·f'(x)|·(1+|fₘₐₓ/f(x)-1|) + 1 ] · ε .
// For x → ∞, we have f(x) → fₘₐₓ which brings us back to
//   |δx/x|  ≈  [ |f(x)|/|x·f'(x)| + 1 ] · ε
// despite the fact that the *relative* accuracy of the complementary value ȳ = f̄(x) diverges according to |δȳ/ȳ|  =  (1+|y/ȳ|)·ε when ȳ → 0 in the limit of x → ∞.
// 
// The attainable *relative* accuracy of 𝛽 = b(s) when s has *relative* accuracy ε is (to lowest order) (|s·b'(s)/b(s)|+1)·ε --- see the source code for a detailed derivation.
// The attainable *relative* accuracy of s = b⁻¹(𝛽) when 𝛽 has *relative* accuracy ε is (to lowest order) (|b(s)/(s·b'(s))|+1)·ε .
// This function returns (s·∂b(x,s)/∂s)/b(x,s,θ=±1). In order to get the accuracy limit of implied volatility calculations, take (1+1/BlackAccuracyFactor(x,s,θ))·DBL_EPSILON.

// NOTE: this function requires s > 0 and θx != 0.
inline double scaled_normalised_black(double θx, double s) {
  assert(s > 0 && θx != 0);
  return (θx > 0 ? normalised_intrinsic(θx) * SQRT_TWO_PI * exp(0.5 * (square(θx / s) + 0.25 * s * s)) : 0) + std::get<0>(scaled_normalised_black_and_ln_vega(-fabs(θx), s));
}

// Returns s · ∂(b(x,s)/∂s) / b(x,s) .
double BlackAccuracyFactor(double x, double s, double θ /* θ=±1 */) {
  // When x = 0, we have bx(x,s) = b(x,s) / (∂(b(x,s)/∂s)  =  s·(1+s²/12+s⁴/240+O(s⁶)) for small s.
  if (0 == x)
    return fabs(s) < DBL_EPSILON ? 1 : s / (erf_cody((0.5 / SQRT_TWO) * s) * SQRT_TWO_PI * exp(0.125 * s * s));
  const double θx = θ < 0 ? -x : x;
  // bx(x,s) = b(x,s) / (∂(b(x, s) / ∂s)
  // For x < 0 (strictly), we have the s→0 asymptotic bx(x,s) = b(x,s) / (∂(b(x, s) / ∂s) → s²/x².
  return s <= 0 ? (θx > 0 ? 0 : DBL_MAX) : s / scaled_normalised_black(θx, s);
}

// Returns  DBL_EPSILON · (1 + |b(x,s) / (s·∂(b(x,s)/∂s)|) .
double ImpliedVolatilityAttainableAccuracy(double x, double s, double θ /* θ=±1 */) {
  // When x = 0, we have bx(x,s) = b(x,s) / (∂(b(x,s)/∂s)  =  s·(1+s²/12+s⁴/240+O(s⁶)) for small s.
  if (0 == x)
    return DBL_EPSILON * (1 + fabs(s <= DBL_EPSILON ? 1 : (erf_cody((0.5 / SQRT_TWO) * s) * SQRT_TWO_PI * exp(0.125 * s * s)) / s));
  const double θx = θ < 0 ? -x : x;
  // For x < 0 (strictly), we have the s→0 asymptotic bx(x,s) = b(x,s) / (∂(b(x, s) / ∂s) → s²/x².
  if (s <= 0) // bx_over_s = (θx > 0 ? DBL_MAX : 0); return std::min( DBL_EPSILON * (1 + bx_over_s) , 1.0 );
    return θx > 0 ? 1 : DBL_EPSILON;
  const double bx = scaled_normalised_black(θx, s);
  return bx * normalised_vega(θx, s) >= DBL_MIN ? DBL_EPSILON * (1 + fabs(bx / s)) : 1;
}
