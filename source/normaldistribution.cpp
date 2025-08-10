//
// This source code resides at www.jaeckel.org/LetsBeRational.7z .
//
// ======================================================================================
//
// normaldistribution.cpp
//

#if defined(_MSC_VER)
# define NOMINMAX // to suppress MSVC's definitions of min() and max()
#endif

#include "normaldistribution.h"
#include <float.h>

namespace {
  // The asymptotic expansion  Φ(z) = φ(z)/|z|·[1-1/z^2+...],  Abramowitz & Stegun (26.2.12), suffices for Φ(z) to have
  // relative accuracy of 1.64E-16 for z<=-10 with 17 terms inside the square brackets (not counting the leading 1).
  // This translates to a maximum of about 9 iterations below, which is competitive with a call to erfc() and never
  // less accurate when z<=-10. Note that, as mentioned in section 4 (and discussion of figures 2 and 3) of George
  // Marsaglia's article "Evaluating the Normal Distribution" (available at http://www.jstatsoft.org/v11/a05/paper),
  // for values of x approaching -8 and below, the error of any cumulative normal function is actually dominated by
  // the hardware (or compiler implementation) accuracy of exp(-x²/2) which is not reliably more than 14 digits when
  // x becomes large. Still, we should switch to the asymptotic only when it is beneficial to do so.
  const double norm_cdf_asymptotic_expansion_first_threshold = -10.0;
  const double norm_cdf_asymptotic_expansion_second_threshold = -1 / sqrt(DBL_EPSILON);
}

double norm_cdf(double z) {
  if (z <= norm_cdf_asymptotic_expansion_first_threshold) {
    // Asymptotic expansion for very negative z following (26.2.12) on page 408
    // in M. Abramowitz and A. Stegun, Pocketbook of Mathematical Functions, ISBN 3-87144818-4.
    double sum = 1;
    if (z >= norm_cdf_asymptotic_expansion_second_threshold) {
      double zsqr = z * z, i = 1, g = 1, x, y, a = DBL_MAX, lasta;
      do {
        lasta = a;
        x = (4 * i - 3) / zsqr;
        y = x * ((4 * i - 1) / zsqr);
        a = g * (x - y);
        sum -= a;
        g *= y;
        ++i;
        a = fabs(a);
      } while (lasta > a && a >= fabs(sum * DBL_EPSILON));
    }
    return -norm_pdf(z) * sum / z;
  }
  return 0.5 * erfc_cody(-z * (1/SQRT_TWO));
}

#if defined( USE_ALGORITHM_AS241 )

//
// ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
//
// Produces the normal deviate Z corresponding to a given lower
// tail area of u; Z is accurate to about 1 part in 10**16.
// see http://lib.stat.cmu.edu/apstat/241 or
// https://csg.sph.umich.edu/abecasis/gas_power_calculator/algorithm-as-241-the-percentage-points-of-the-normal-distribution.pdf.
//

namespace {

  const double split1 = 0.425;
  const double split2 = 5.0;
  const double const1 = 0.180625; // = split1²
  const double const2 = 1.6;

  // Coefficients for P close to 0.5
  const double A0 = 3.3871328727963666080E0;
  const double A1 = 1.3314166789178437745E+2;
  const double A2 = 1.9715909503065514427E+3;
  const double A3 = 1.3731693765509461125E+4;
  const double A4 = 4.5921953931549871457E+4;
  const double A5 = 6.7265770927008700853E+4;
  const double A6 = 3.3430575583588128105E+4;
  const double A7 = 2.5090809287301226727E+3;
  const double B1 = 4.2313330701600911252E+1;
  const double B2 = 6.8718700749205790830E+2;
  const double B3 = 5.3941960214247511077E+3;
  const double B4 = 2.1213794301586595867E+4;
  const double B5 = 3.9307895800092710610E+4;
  const double B6 = 2.8729085735721942674E+4;
  const double B7 = 5.2264952788528545610E+3;
  // Coefficients for P not close to 0, 0.5 or 1.
  const double C0 = 1.42343711074968357734E0;
  const double C1 = 4.63033784615654529590E0;
  const double C2 = 5.76949722146069140550E0;
  const double C3 = 3.64784832476320460504E0;
  const double C4 = 1.27045825245236838258E0;
  const double C5 = 2.41780725177450611770E-1;
  const double C6 = 2.27238449892691845833E-2;
  const double C7 = 7.74545014278341407640E-4;
  const double D1 = 2.05319162663775882187E0;
  const double D2 = 1.67638483018380384940E0;
  const double D3 = 6.89767334985100004550E-1;
  const double D4 = 1.48103976427480074590E-1;
  const double D5 = 1.51986665636164571966E-2;
  const double D6 = 5.47593808499534494600E-4;
  const double D7 = 1.05075007164441684324E-9;
  // Coefficients for P very close to 0 or 1
  const double E0 = 6.65790464350110377720E0;
  const double E1 = 5.46378491116411436990E0;
  const double E2 = 1.78482653991729133580E0;
  const double E3 = 2.96560571828504891230E-1;
  const double E4 = 2.65321895265761230930E-2;
  const double E5 = 1.24266094738807843860E-3;
  const double E6 = 2.71155556874348757815E-5;
  const double E7 = 2.01033439929228813265E-7;
  const double F1 = 5.99832206555887937690E-1;
  const double F2 = 1.36929880922735805310E-1;
  const double F3 = 1.48753612908506148525E-2;
  const double F4 = 7.86869131145613259100E-4;
  const double F5 = 1.84631831751005468180E-5;
  const double F6 = 1.42151175831644588870E-7;
  const double F7 = 2.04426310338993978564E-15;
  
#define AB(r) (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0) / (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + 1.0)
#define CD(r) (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0) / (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + 1.0)
#define EF(r) (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0) / (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + 1.0)

}

double inverse_norm_cdf(double u) {
#ifdef THIS_IS_NOT_NEEDED // The result will be NaN if u<=0 or u>=1 with or without these early checks.
  if (u <= 0)
    return log(u);
  if (u >= 1)
    return log(1 - u);
#endif
  const double q = u - 0.5;
  if (fabs(q) <= split1) {
    const double r = const1 - q * q;
    return q * AB(r);
  } else {
    double r = q < 0.0 ? u : 1.0 - u;
    r = sqrt(-log(r));
    double ret;
    if (r < split2) {
      r = r - const2;
      ret = CD(r);
    } else {
      r = r - split2;
      ret = EF(r);
    }
    return q < 0.0 ? -ret : ret;
  }
}

// We can use the internal branches of Φ⁻¹(·) to implement erfinv() avoiding catastrophic subtractive cancellation for small arguments.
double erfinv(double 𝑒) {
  // Φ(x) = erfc(-x/√2)/2 = (1+erf(x/√2))/2 = erf(x/√2)/2 + half.
  // erf(z) = 2 · ( Φ(√2·z) - half ).
  // Hence, if 𝑒 = erfc(z), and y is the solution to Φ(y) - half = 𝑒/2, then z = y/√2.
  const double q = 0.5 * 𝑒;
  if (fabs(q) <= split1) {
    const double r = const1 - q * q;
    return q * AB(r) * (1 / SQRT_TWO);
  }
  double r = -0.5 * fabs(𝑒) + 0.5;
  r = sqrt(-log(r));
  double ret;
  if (r < split2) {
    r = r - const2;
    ret = CD(r);
  } else {
    r = r - split2;
    ret = EF(r);
  }
  return (q < 0.0 ? -ret : ret) * (1 / SQRT_TWO);
}

#else

#include <assert.h>

//
// This source code resides at www.jaeckel.org/LetsBeRational.7z .
//
// ======================================================================================
// Copyright © 2024 Peter Jäckel.
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                     //
// PLEASE REFER TO THE BELOW ALGORITHM FOR INVERSE NORMAL PROBABILITY AS "PJ-2024-Inverse-Normal"      //
//                                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////////////////////


#if !defined( USE_SPECIALISED_RATIONAL_FUNCTION_APPROXIMATION_FOR_ATM_IMPLIED_VOLATILITY )
inline
#endif
// Specialisation of x = Φ⁻¹(p) for x ≤ -1, p = Φ(x) ≤ 0.1586552539314570514148
double inverse_norm_cdf_for_low_probabilities(double p) {
  assert(!(p > 0.15865525393146)); // This formulation permits NaN to pass through below. The result will be NaN as it should.
  // Five branches.
  //                    │          I.       │       II.       │      III.      │      IV.      │        V.          │
  //                    │         ←-→       │       ←-→       │      ←-→       │      ←-→      │       ←-→          │
  //  In r=√|ln(p)|:    [ 1.35684  ,      2.05       ,      3.41      ,       6.7      ,     12.9       ,    27.314 ]
  //  In |x|:           [    1     ,     2.1712      ,     4.2905     ,     9.1375     ,    18.0331     ,    38.509 ]
  //  In p = Φ(-|x|):   [ 0.15866  ,    0.014958     ,    8.912E-6    ,    3.1954E-20  ,   5.35865E-73  ,    1E-324 ]   ==>  Compare this with DBL_TRUE_MIN = 4.94E-324.
  const double r = sqrt(-log(p));
  // All of the below are Remez-optimized minimax rational function approximations of order (5,5).
  if (r < 6.7)
    if (r < 3.41)
      if (r < 2.05) // I.   Accuracy better than 7.6E-17 in perfect arithmetic.
        return (3.691562302945566191 + r * (4.7170590600740689449E1 + r * (6.5451292110261454609E1 + r * (-7.4594687726045926821E1 + r * (-8.3383894003636969722E1 - 1.3054072340494093704E1 * r))))) / (1 + r * (2.0837211328697753726E1 + r * (7.1813812182579255459E1 + r * (5.9270122556046077717E1 + r * (9.2216887978737432303 + 1.8295174852053530579E-4 * r)))));
      else          // II.  Accuracy better than 9.4E-17 in perfect arithmetic.
        return (3.2340179116317970288 + r * (1.449177828689122096E1 + r * (6.8397370256591532878E-1 + r * (-1.81254427791789183E1 + r * (-1.005916339568646151E1 - 1.2013147879435525574E0 * r))))) / (1 + r * (8.8820931773304337525 + r * (1.4656370665176799712E1 + r * (7.1369811056109768745 + r * (8.4884892199149255469E-1 + 1.0957576098829595323E-5 * r)))));
    else            // III. Accuracy better than 9.1E-17 in perfect arithmetic.
      return (3.1252235780087584807 + r * (9.9483724317036560676 + r * (-5.1633929115525534628 + r * (-1.1070534689309368061E1 + r * (-2.8699061335882526744 - 1.5414319494013597492E-1 * r))))) / (1 + r * (7.076769154309171622 + r * (8.1086341122361532407 + r * (2.0307076064309043613 + r * (1.0897972234131828901E-1 + 1.3565983564441297634E-7 * r)))));
  else
    if (r < 12.9)   // IV.  Accuracy better than 9E-17 in perfect arithmetic.
      return (2.6161264950897283681 + r * (2.250881388987032271 + r * (-3.688196041019692267 + r * (-2.9644251353150605663 + r * (-4.7595169546783216436E-1 - 1.612303318390145052E-2 * r))))) / (1 + r * (3.2517455169035921495 + r * (2.1282030272153188194 + r * (3.3663746405626400164E-1 + r * (1.1400087282177594359E-2 + 3.0848093570966787291E-9 * r)))));
    else            // V.   Accuracy better than 9.5E-17 in perfect arithmetic
      return (2.3226849047872302955 + r * (-4.2799650734502094297E-2 + r * (-2.5894451568465728432 + r * (-8.6385181219213758847E-1 + r * (-6.5127593753781672404E-2 - 1.0566357727202585402E-3 * r))))) / (1 + r * (1.9361316119254412206 + r * (6.1320841329197493341E-1 + r * (4.6054974512474443189E-2 + r * (7.471447992167225483E-4 + 2.3135343206304887818E-11 * r)))));
}

namespace {
  // uₘₐₓ = Φ(1) - half = 0.3413447460685429
  const double uₘₐₓ = 0.3413447460685429;
}

// Inverse of Φ(x)-half, i.e., for a given u ∈ [-uₘₐₓ,uₘₐₓ] find x such that u = Φ(x)-half.
inline double inverse_norm_cdfmhalf_for_midrange_probabilities(double u) {
  assert(!(fabs(u) > uₘₐₓ)); // This formulation permits NaN to pass through below. The result will be NaN as it should.
  // innermost branch: x ∈ [-1,1], since u = Φ(x)-half ⟹ u ∈ [-uₘₐₓ,uₘₐₓ] with uₘₐₓ=Φ(1)-half=0.3413447460685429, and p ∈ [0.15866,0.84134].
    // Define f(x) := √(2π)·(Φ(x)-half). It is an odd function: f(-x) = -f(x). Its Taylor expansion is f(x) = x - x³/6 + x⁵/40 - O(x⁷).
    // Set y = f(x). The inverse function f⁻¹(y) is also odd with its Taylor expansion as f⁻¹(y) = y + y³/6 + 7·y⁵/120 + O(y⁷).
    // Define g̃(y) := (f⁻¹(y)/y-1)/y² which is now an even function and we can focus on y ≥ 0. We define q := y² and g(q) := g̃(√q) = (f⁻¹(√q)/√q-1)/q.
    // We then compute a Remez-optimized minimax rational function approximation of order (5,5) for g(q).
    // In the end, we transform to h(s) := √(2π)·(1+2π·(uₘₐₓ²-s)·g(2π·(uₘₐₓ²-s)))   with   s = uₘₐₓ² - u²  where uₘₐₓ = f(xₘₐₓ)/√(2π) = (Φ(xₘₐₓ)-half).  Note that this makes h(s) a (6,5) rational function.
    // Accuracy better than 9.8E-17 in perfect arithmetic within this branch.
  const double s = uₘₐₓ * uₘₐₓ - u * u;
  return u * ((2.92958954698308805 + s * (5.0260572167303103E1 + s * (3.01870541922933937E2 + s * (7.4997781456657924E2 + s * (6.90489242061408612E2 + s * (1.34233243502653864E2 - 7.58939881401259242 * s)))))) / (1 + s * (1.8918538074574598E1 + s * (1.29404120448755281E2 + s * (3.86821208540417453E2 + s * (4.79123914509756757E2 + 1.79227008508102628E2 * s))))));
}

// Inverse of Φ(x). Given p, return x such that p = Φ(x), i.e., x = Φ⁻¹(p).
double inverse_norm_cdf(double p) {
  const double u = p - 0.5;
  if (fabs(u) < uₘₐₓ)
      return inverse_norm_cdfmhalf_for_midrange_probabilities(u);
  // Tail probability min(p,1-p) = 0.15866. r = √(-ln(min(p,1-p))) > 1.3568425277.
  return u > 0 ? -inverse_norm_cdf_for_low_probabilities(1 - p) : inverse_norm_cdf_for_low_probabilities(p);
}

// We can use the internal branches of Φ⁻¹(·) to implement erfinv() avoiding catastrophic subtractive cancellation for small arguments.
double erfinv(double 𝑒) {
  // Φ(x) = erfc(-x/√2)/2 = (1+erf(x/√2))/2 = erf(x/√2)/2 + half.
  // erf(z) = 2 · ( Φ(√2·z) - half ).
  // Hence, if 𝑒 = erfc(z), and y is the solution to Φ(y) - half = 𝑒/2, then z = y/√2.
  if (fabs(𝑒) < (2 * uₘₐₓ))
    return inverse_norm_cdfmhalf_for_midrange_probabilities(0.5 * 𝑒) * (1 / SQRT_TWO);
  return (𝑒 < 0 ? inverse_norm_cdf_for_low_probabilities(0.5 * 𝑒 + 0.5) : -inverse_norm_cdf_for_low_probabilities(-0.5 * 𝑒 + 0.5)) * (1 / SQRT_TWO);
}

#endif
