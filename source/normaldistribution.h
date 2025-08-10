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
#ifndef   NORMAL_DISTRIBUTION_H
#define   NORMAL_DISTRIBUTION_H

#include <math.h>
#include <cmath>

#define SQRT_TWO              1.4142135623730950488016887242096980785696718753769
#define SQRT_TWO_PI           2.5066282746310005024157652848110452530069867406099
#define LN_TWO_PI             1.8378770664093454835606594728112352797227949472756

double erf_cody(double z);
double erfc_cody(double z);
double erfcx_cody(double z);
double erfinv(double 𝑒);

inline double norm_pdf(double x) { return (1 / SQRT_TWO_PI) * exp(-.5 * x * x); }

extern "C" double norm_cdf(double z);
extern "C" double inverse_norm_cdf(double u);

#endif // NORMAL_DISTRIBUTION_H
