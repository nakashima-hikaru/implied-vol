//
// This source code resides at www.jaeckel.org/ImpliedNormalVolatility.7z .
//
// ======================================================================================
// Copyright © 2022 Peter Jäckel.
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

/////////////////////////////////////////////////////////////////////////////////////////
//
// This is a reference implementation of the analytical normal (Bachelier) implied
// volatility function by Peter Jaeckel, © 2016-2022.
//
// See www.jaeckel.org/ImpliedNormalVolatility.pdf for details of the mathematics,
// published January 2017; Wilmott, pages 52-54, March 2017.
//
/////////////////////////////////////////////////////////////////////////////////////////

#define NOMINMAX // to suppress MSVC's definition of min() and max()
#include <float.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <assert.h>

#include "normaldistribution.h"

constexpr double ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;

// See http://gcc.gnu.org/wiki/Visibility
#ifndef EXPORT_EXTERN_C
# if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
#  define EXPORT_EXTERN_C extern "C" __declspec(dllexport)
# else
#  define EXPORT_EXTERN_C extern "C" __attribute__((visibility("default")))
# endif
#endif

namespace {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
}

EXPORT_EXTERN_C double NormPdf(double x) { return norm_pdf(x); }
EXPORT_EXTERN_C double NormCdf(double x) { return norm_cdf(x); }

// We use an implementation of (Φ̃(x)·x) = Φ(x)·x+φ(x) for the [somewhat indirect] evaluation of Φ̃(x),
// but a direct implementation of the inverse of Φ̃(x) for numerical efficiency reasons.
double PhiTildeTimesX(double x);

//
// Φ̃(x) := Φ(x)+φ(x)/x --- careful: this function incurs massive subtractive cancellation errors even for just moderately negative x values.
//
// Note that Φ̃(x) has an unusual pole at x=0. Specifically, we have
//
//     Φ̃(x) < 0  when  x < 0     with    Φ̃(x) ->  0  [from below]  when  x -> -∞
//     Φ̃(x) > 1  when  x > 0     with    Φ̃(x) ->  1  [from above]  when  x -> +∞
//
// i.e, different horizontal asymptotic levels on the left and the right,
//
//                  │•
//                  │•          • Φ̃(x)
//                  │•
//                  │ •
//                  │ •
//                  │  ••
//                  1    ••••••••••••••
//                  │
//                  │
//                  │
//  ••••••••••••─── 0 ──────────────────> x
//              ••  │
//                • │
//                • │
//                 •│
//                 •│
//                 •│
//             
// unlike, e.g., 1/x or similar. The upper right quadrant of Φ̃(x) is transformed into the lower left via
//
//     Φ̃(x) = 1-Φ̃(-x).
//
// For large negative x, this function incurs two fundamental numerical problems:
//
//  1) It converges so rapidly to zero that for x < -37.326030040668367743 we have |Φ̃(x)| < DBL_MIN
//     with DBL_MIN = 2.2250738585072014E-308.
//
//  1) It incurs a subtractive cancellation error as can be seen from Abramowitz-Stegun (26.2.12) for x → -∞
//
//         Φ(x)   ≈   -φ(x)/x·( 1 - 1/x^2 + 3/x^4 - 3·5/x^6 + 3·5·7/x^6 - ... ) .
//
//     In fact, since
//                       Φ̃(x) = (  1  +  x·Φ(x)/φ(x)  )  ·  φ(x)/x
//     and
//              x·Φ(x)/φ(x) < -1/2    when    x ≤ -0.61200318096248076056 =: x_c
//
//     we have subtractive cancellation as soon as x < x_c as defined above.
//
// Notes:
//
//  +) Φ̃(x_c) = -0.27026782648771585618
//
//  +) ln|Φ̃(x_c)| = -1.3083418616865091797
//
//  +) 0 > Φ̃(x) > -1  for  x < -0.2760298047981433 and thus sqrt(-ln(-Φ̃(x))) is well defined when x < -0.2760298047981433
//
EXPORT_EXTERN_C double PhiTilde(double x) {
#ifdef USE_SIMPLISTIC_VERSION_OF_PHI_TILDE
  // Equation (1.2)
  return norm_cdf(x)+norm_pdf(x)/x;
#else
  return PhiTildeTimesX(x) / x;
#endif
}

//
// f(x) := Φ(x)·x+φ(x).
//
// Note that ƒ(x) := f(x)-x/2 is a symmmettric function, i.e., ƒ(x)=ƒ(-x), which permits expansion in x².
//
// Also note that both Φ(x)·x+φ(x) and Φ̃(x):=Φ(x)+φ(x)/x incur loss of precision due to subtractive cancellation
// as and when Φ(x)·x < 0 such that Φ(x)·x < -φ(x)/2. This happens when Φ(x_c)·x_c/φ(x_c) = -1/2, the solution to
// which is x_c = -0.61200318096248076056.
//
//
//  Note that f(x) = f(-x)+x which we use for x > |x_c|.
//
double PhiTildeTimesX(double x) {
#ifdef USE_SIMPLISTIC_VERSION_OF_PHI_TILDE_TIMES_X
  return norm_cdf(x)*x+norm_pdf(x);
#else
  // We follow the footsteps of Cody in his implementation of erf(), erfc() and erfcx() as published in "Rational
  // Chebyshev approximations for the error function", W. J. Cody, Math. Comp., 1969, pp. 631-638, by dividing
  // into three (slightly differently sized) branches, and using (essentially) the same transformations and symmetries.
  // In the first branch, we directly use a linear Chebyshev-Padé approximant. In the second and third branch, the
  // respective approximant was improved with a Remez-II implementation to arrive at full (standard IEEE 754
  // with 53 bit mantissa) double floating point precision.
  if (fabs(x) <= 0.61200318096248076056) {
    // -0.61200318096248076056 ≤ x ≤ 0.61200318096248076056
    // Linear Chebyshev-Padé approximant for g := ( φ(x) + x·Φ(x) - φ(0) - x·Φ(0) ) / x² = ( 1 - w/12 + w²/120 - w³/1344 + ... ) / √(8π) with w := x².
    const double h = (x * x - 1.8727394675409748661E-1) * 5.3397710537550806412;
    const double g = (1.9641549843774702457E-1 + h * (2.9444812226268915305E-3 + 3.095828855856470717E-5 * h)) / (1 + h * (3.0261016846592326803E-2 + h * (3.3735461911896198861E-4 + h * (1.290112376540573289E-6 - 1.6711975835244204502E-9 * h))));
    // |∆f/f-1| < 2.5E-18
    return ONE_OVER_SQRT_TWO_PI + x*(0.5+x*g);
  }
  if (x > 0)
    return PhiTildeTimesX(-x)+x;
  if (x >= -3.5) {
    // -3.5 ≤ x < -0.61200318096248076056
    // Remez-II optimized minimax rational function based on linear Chebyshev-Padé approximant initial guess for g(x) := exp(x²/2)·f(x).
    // |∆g/g-1| < 1.6E-20,     |∆f/f-1| = |∆g/g-1| (in perfect arithmetic when the below exponential incurs no error).
    const double g = (3.9894228040096173296E-1 + x * ((-2.8827250122716400843E-1) + x * (1.1748934770055073669E-1 + x * ((-2.9208930498324232842E-2) + x * (4.6704817087348921557E-3 + x * ((-4.4448405482476358857E-4) + x * (1.9865267442385935787E-5 + x * (7.6387393474143610035E-10 + 1.3291525220137582449E-11 * x)))))))) / (1 + x * ((-1.9759061396728604494) + x * (1.7709332198933623888 + x * ((-9.4350250026446231963E-1) + x * (3.2816118145388593816E-1 + x * ((-7.6697408088214742324E-2) + x * (1.1843224303096222834E-2 + x * ((-1.1151416365524860908E-3) + 4.9741005333758689307E-5 * x))))))));
    return exp(-0.5*(x*x))*g;
  }
  // x < -3.5
  // Linear Chebyshev-Padé approximant for g :=  x²·(1-x²·f(x)/φ(x)) as a function of w := 1/x².
  double w = 1 / (x * x); // w < 0.081632653
  const double g = (2.999999999999991221 + w * (2.3654556627823149931E2 + w * (6.8126773449358787324E3 + w * (8.9697941598360784061E4 + w * (5.5163920591268613879E5 + w * (1.4345061123335662019E6 + w * (1.1504988246344881836E6 + 1.1867600400997691371E4 * w))))))) / (1 + w * (8.3848522092737134602E1 + w * (2.6551350587809577877E3 + w * (4.0555290884673789153E4 + w * (3.166737476299376429E5 + w * (1.2329795958024320559E6 + w * (2.1409810540619049948E6 + 1.2145667804093160403E6 * w)))))));
  // f(x) = φ(x)/x²·(1-g/x²) = φ(x)·w·(1-g·w) . Note that g(w) is a number between two and three in this interval, and g(w=0) = 3.
  // |∆f/f-1| < 2.5E-18 in -38.5 ≤ x ≤ -3.5, which corresponds to 6.74650025299375948726598E-4 ≤ w ≤ 8.163265306122448979591836b-2.
  // Note that |f(x)| < φ(x)/x² and that φ(x)/x² < (DBL_MIN·DBL_EPSILON), i.e., underflows to zero, when x < -38.37250105526059780836052 .
  return ONE_OVER_SQRT_TWO_PI * exp(-0.5*(x*x)) * w * (1 - g * w);
#endif
}

EXPORT_EXTERN_C double InvPhiTilde(double phi_tilde_star) {
  if (phi_tilde_star > 1)
    return -InvPhiTilde(1 - phi_tilde_star);
  if (phi_tilde_star >= 0)
    return NaN;
  double x_bar = NaN;
  if (phi_tilde_star < -0.00188203927) {
    // Equation (2.1)
    const double g = 1/(phi_tilde_star-0.5), g2 = g*g;
    // Equation (2.2)
    const double xi_bar = (0.032114372355-g2*(0.016969777977-g2*(0.002620733246-0.000096066952861*g2)))/(1-g2*(0.6635646938-g2*(0.14528712196-0.010472855461*g2)));
    // Equation (2.3)
    x_bar = g*(ONE_OVER_SQRT_TWO_PI+xi_bar*g2);
  } else {
    // Equation (2.4)
    const double h = sqrt(-log(-phi_tilde_star));
    // Equation (2.5)
    x_bar = (9.4883409779-h*(9.6320903635-h*(0.58556997323+2.1464093351*h)))/(1-h*(0.65174820867+h*(1.5120247828+0.000066437847132*h)));
  }
  // Equation (2.7)
  const double q = (PhiTilde(x_bar)-phi_tilde_star)/norm_pdf(x_bar), x2 = x_bar*x_bar;
  // Equation (2.6)
  const double x_star = x_bar + 3*q*x2*(2-q*x_bar*(2+x2))/(6+q*x_bar*(-12+x_bar*(6*q+x_bar*(-6+q*x_bar*(3+x2)))));
  return x_star;
  
}

inline double intrinsic_value(double forward, double strike, double q /* q=±1 */){
  return fabs(std::max((q<0?strike-forward:forward-strike),0.0));
}

EXPORT_EXTERN_C double Bachelier(double forward, double strike, double sigma, double T, double q /* q=±1 */) {
  double s = fabs(sigma) * sqrt(T);
  if (s<DBL_MIN)
    return intrinsic_value(forward, strike, q);
  double theta = q<0?-1:1, moneyness = theta*(forward-strike), x = moneyness/s;
  return s*PhiTildeTimesX(x);
}

EXPORT_EXTERN_C double ImpliedNormalVolatility(double price, double forward, double strike, double T, double q /* q=±1 */) {
  if (forward == strike)
    return price * SQRT_TWO_PI / sqrt(T);
  const double intrinsic = intrinsic_value(forward, strike, q), absolute_moneyness = fabs(forward - strike);
  if (price == intrinsic)
    return 0.0;
  if (price < intrinsic)
    // This signals that the price is below the intrinsic value.
    return -DBL_MAX;
  // Equation (1.6)
  const double phi_tilde_star = (intrinsic - price) / absolute_moneyness;
  // Solve equation (1.7)
  const double x_star = InvPhiTilde(phi_tilde_star);
  // Equation (1.8)
  return absolute_moneyness / fabs(x_star * sqrt(T));
}

// NOTE: any (worksheet) function name specified with a leading '.' will be registered preceded by this XLL's name.

// An alternative syntax to the above, which needs C++-11, would be:
//  const char *ImpliedNormalVolatilitySynonyms[] = { ".ImpliedBachelierVolatility", ".ImpliedNormalVolatility", "ImpliedBachelierVolatility", "ImpliedNormalVolatility", 0 };
//  DECLARE_XL_FUNCTION_WITH_SYNOYMS( ImpliedNormalVolatility, ImpliedNormalVolatilitySynonyms, "BBBBBB", "price,forward,strike,T,q", "returns the normal (aka 'Bachelier') volatility that matches the given price.", "the price", "the forward", "the strike", "the time to expiry", "q=±1 for calls/puts")
