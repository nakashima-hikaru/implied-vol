// Error-function primitives adapted from Peter Jäckel's Let's Be Rational 1.0.0.1520.
// Copyright © 2013-2024 Peter Jäckel.
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
// The Software is provided "as is" without warranty of any kind, either express or
// implied, including without limitation any implied warranties of condition,
// uninterrupted use, merchantability, fitness for a particular purpose, or non-infringement.
// Cody's original error-function formulation: http://www.netlib.org/specfun/erf
// Preserve the archived decimal constants and their original rounding.
#![allow(clippy::excessive_precision, clippy::approx_constant)]

#[inline]
fn ab(z: f64) -> f64 {
    const A: [f64; 5] = [
        3.1611237438705656,
        113.864154151050156,
        377.485237685302021,
        3209.37758913846947,
        0.185777706184603153,
    ];
    const B: [f64; 4] = [
        23.6012909523441209,
        244.024637934444173,
        1282.61652607737228,
        2844.23683343917062,
    ];
    ((((A[4] * z + A[0]) * z + A[1]) * z + A[2]) * z + A[3])
        / ((((z + B[0]) * z + B[1]) * z + B[2]) * z + B[3])
}

#[inline]
fn cd(y: f64) -> f64 {
    const C: [f64; 9] = [
        0.564188496988670089,
        8.88314979438837594,
        66.1191906371416295,
        298.635138197400131,
        881.95222124176909,
        1712.04761263407058,
        2051.07837782607147,
        1230.33935479799725,
        2.15311535474403846e-8,
    ];
    const D: [f64; 8] = [
        15.7449261107098347,
        117.693950891312499,
        537.181101862009858,
        1621.38957456669019,
        3290.79923573345963,
        4362.61909014324716,
        3439.36767414372164,
        1230.33935480374942,
    ];
    let mut num = C[8] * y + C[0];
    let mut den = y + D[0];
    for i in 1..8 {
        num = num * y + C[i];
        den = den * y + D[i];
    }
    num / den
}

#[inline]
fn pq(z: f64) -> f64 {
    const P: [f64; 6] = [
        0.305326634961232344,
        0.360344899949804439,
        0.125781726111229246,
        0.0160837851487422766,
        6.58749161529837803e-4,
        0.0163153871373020978,
    ];
    const Q: [f64; 5] = [
        2.56852019228982242,
        1.87295284992346047,
        0.527905102951428412,
        0.0605183413124413191,
        0.00233520497626869185,
    ];
    z * (((((P[5] * z + P[0]) * z + P[1]) * z + P[2]) * z + P[3]) * z + P[4])
        / (((((z + Q[0]) * z + Q[1]) * z + Q[2]) * z + Q[3]) * z + Q[4])
}

#[inline]
pub(crate) fn erfcx_cody(x: f64) -> f64 {
    let y = x.abs();
    if y <= 0.46875 {
        let z = y * y;
        return z.exp() * (1.0 - x * ab(z));
    }
    if x < -26.6287357137514 {
        return f64::MAX;
    }
    let result = if y <= 4.0 {
        cd(y)
    } else {
        (0.56418958354775628695 - pq(1.0 / (y * y))) / y
    };
    if x < 0.0 {
        let xt = (x * 16.0).trunc() / 16.0;
        let expx2 = (xt * xt).exp() * ((x - xt) * (x + xt)).exp();
        return (expx2 + expx2) - result;
    }
    result
}

#[inline]
pub(crate) fn erfc_cody(x: f64) -> f64 {
    let y = x.abs();
    if y <= 0.46875 {
        return 1.0 - x * ab(y * y);
    }
    let erfc_abs = if y >= 26.543 {
        0.0
    } else {
        let scaled = if y <= 4.0 {
            cd(y)
        } else {
            (0.56418958354775628695 - pq(1.0 / (y * y))) / y
        };
        let yt = (y * 16.0).trunc() / 16.0;
        scaled * ((-yt * yt).exp() * (-(y - yt) * (y + yt)).exp())
    };
    if x < 0.0 { 2.0 - erfc_abs } else { erfc_abs }
}

#[inline(always)]
fn inverse_stage<const FUSED: bool>(a: f64, b: f64, c: f64) -> f64 {
    if FUSED { a.mul_add(b, c) } else { a * b + c }
}

#[inline]
pub(crate) fn erfinv(e: f64) -> f64 {
    erfinv_impl::<true>(e)
}

#[inline]
pub(crate) fn erfinv_unfused(e: f64) -> f64 {
    erfinv_impl::<false>(e)
}

#[inline]
pub(crate) fn inverse_norm_cdf(p: f64) -> f64 {
    inverse_norm_cdf_impl::<true>(p)
}

#[inline]
pub(crate) fn inverse_norm_cdf_unfused(p: f64) -> f64 {
    inverse_norm_cdf_impl::<false>(p)
}

#[inline]
fn erfinv_impl<const FUSED: bool>(e: f64) -> f64 {
    // Multiplication order follows the archived PJ-2024 implementation.
    const INV_SQRT_TWO: f64 = 1.0 / 1.4142135623730950488016887242096980785696718753769;
    if e.abs() < 2.0 * 0.3413447460685429 {
        return inverse_norm_mid::<FUSED>(0.5 * e) * INV_SQRT_TWO;
    }
    (if e < 0.0 {
        inverse_norm_low::<FUSED>(0.5 * e + 0.5)
    } else {
        -inverse_norm_low::<FUSED>(-0.5 * e + 0.5)
    }) * INV_SQRT_TWO
}

// Share coefficients and branches while selecting Horner rounding at compile time.
// Main solver seeds use FMA. The original LBR fallback relies on its unfused seed
// rounding: fusing it changed later LBR iterates and worsened some small-a roots.
// Both choices are monomorphized; there is no runtime mode branch.
// Copyright © 2024 Peter Jäckel.
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
// WARRANTY DISCLAIMER
// The Software is provided "as is" without warranty of any kind, either express or implied,
// including without limitation any implied warranties of condition, uninterrupted use,
// merchantability, fitness for a particular purpose, or non-infringement.
// Please refer to these rational expressions as PJ-2024-Inverse-Normal.
#[inline]
fn pj2024_tail_from_r<const FUSED: bool>(r: f64) -> f64 {
    if r < 2.05 {
        inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            -f64::from_bits(0x402a1baf5eac14a3),
                            r,
                            -f64::from_bits(0x4054d891b827b015),
                        ),
                        -f64::from_bits(0x4052a60f5d1baccb),
                    ),
                    f64::from_bits(0x40505ce1f84da102),
                ),
                f64::from_bits(0x404795d5e9ad97d7),
            ),
            f64::from_bits(0x400d8851d1126197),
        ) / inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            f64::from_bits(0x3f27fad78da3fba2),
                            r,
                            f64::from_bits(0x4022718131b183ba),
                        ),
                        f64::from_bits(0x404da293603c109e),
                    ),
                    f64::from_bits(0x4051f4157fb150ea),
                ),
                f64::from_bits(0x4034d6537b4c98fa),
            ),
            f64::from_bits(0x3ff0000000000000),
        )
    } else if r < 3.41 {
        inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            -f64::from_bits(0x3ff33895dae6b322),
                            r,
                            -f64::from_bits(0x40241e4aaa232ff4),
                        ),
                        -f64::from_bits(0x4032201d049a17c2),
                    ),
                    f64::from_bits(0x3fe5e31cd17b0745),
                ),
                f64::from_bits(0x402cfbca5d162951),
            ),
            f64::from_bits(0x4009df44c8691823),
        ) / inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            f64::from_bits(0x3ee6facdcaa72548),
                            r,
                            f64::from_bits(0x3feb29c536e6589a),
                        ),
                        f64::from_bits(0x401c8c44c6631387),
                    ),
                    f64::from_bits(0x402d500fd0d9f9eb),
                ),
                f64::from_bits(0x4021c3a1b7895161),
            ),
            f64::from_bits(0x3ff0000000000000),
        )
    } else if r < 6.7 {
        inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            -f64::from_bits(0x3fc3baf6d6959c71),
                            r,
                            -f64::from_bits(0x4006f59158d2cf92),
                        ),
                        -f64::from_bits(0x4026241d1f6fa26f),
                    ),
                    -f64::from_bits(0x4014a75078ae10d5),
                ),
                f64::from_bits(0x4023e591124530be),
            ),
            f64::from_bits(0x400900753821e2c2),
        ) / inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            f64::from_bits(0x3e82353c88a0efe5),
                            r,
                            f64::from_bits(0x3fbbe61857621453),
                        ),
                        f64::from_bits(0x40003ee3a12adf96),
                    ),
                    f64::from_bits(0x4020379ee3ee918c),
                ),
                f64::from_bits(0x401c4e9c92bc65db),
            ),
            f64::from_bits(0x3ff0000000000000),
        )
    } else if r < 12.9 {
        inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            -f64::from_bits(0x3f90828e70f2df2e),
                            r,
                            -f64::from_bits(0x3fde75fe19a0a553),
                        ),
                        -f64::from_bits(0x4007b724867cf2cf),
                    ),
                    -f64::from_bits(0x400d816ced0b547c),
                ),
                f64::from_bits(0x400201ce1a06feaf),
            ),
            f64::from_bits(0x4004edd3ba54e03b),
        ) / inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            f64::from_bits(0x3e2a7f9148ad1807),
                            r,
                            f64::from_bits(0x3f8758edd0635312),
                        ),
                        f64::from_bits(0x3fd58b77dcaeb8cd),
                    ),
                    f64::from_bits(0x4001068f4f091a9c),
                ),
                f64::from_bits(0x400a039327501fd5),
            ),
            f64::from_bits(0x3ff0000000000000),
        )
    } else {
        inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            -f64::from_bits(0x3f514fda059b852b),
                            r,
                            -f64::from_bits(0x3fb0ac33b53d5ac3),
                        ),
                        -f64::from_bits(0x3feba4ac8e3e93fd),
                    ),
                    -f64::from_bits(0x4004b72f05bb8871),
                ),
                -f64::from_bits(0x3fa5e9d5f85eac60),
            ),
            f64::from_bits(0x400294dbd2c7cacf),
        ) / inverse_stage::<FUSED>(
            r,
            inverse_stage::<FUSED>(
                r,
                inverse_stage::<FUSED>(
                    r,
                    inverse_stage::<FUSED>(
                        r,
                        inverse_stage::<FUSED>(
                            f64::from_bits(0x3db970052b2f23e1),
                            r,
                            f64::from_bits(0x3f487b813d2f81e3),
                        ),
                        f64::from_bits(0x3fa7948482b2c933),
                    ),
                    f64::from_bits(0x3fe39f674017131f),
                ),
                f64::from_bits(0x3ffefa65241f8b67),
            ),
            f64::from_bits(0x3ff0000000000000),
        )
    }
}

#[inline]
fn inverse_norm_mid<const FUSED: bool>(u: f64) -> f64 {
    const UMAX: f64 = 0.3413447460685429;
    let s = UMAX * UMAX - u * u;
    u * (inverse_stage::<FUSED>(
        s,
        inverse_stage::<FUSED>(
            s,
            inverse_stage::<FUSED>(
                s,
                inverse_stage::<FUSED>(
                    s,
                    inverse_stage::<FUSED>(
                        s,
                        inverse_stage::<FUSED>(
                            -f64::from_bits(0x401e5b8b5cd9f0eb),
                            s,
                            f64::from_bits(0x4060c776bb13fce4),
                        ),
                        f64::from_bits(0x408593e9f7bdec9f),
                    ),
                    f64::from_bits(0x40876fd290718815),
                ),
                f64::from_bits(0x4072ddedbd5e0cc5),
            ),
            f64::from_bits(0x4049215a6dc46847),
        ),
        f64::from_bits(0x40076fcca4f7f770),
    ) / inverse_stage::<FUSED>(
        s,
        inverse_stage::<FUSED>(
            s,
            inverse_stage::<FUSED>(
                s,
                inverse_stage::<FUSED>(
                    s,
                    inverse_stage::<FUSED>(
                        f64::from_bits(0x40666743a758c6dd),
                        s,
                        f64::from_bits(0x407df1fb8dc7ee7b),
                    ),
                    f64::from_bits(0x40782d23ab91049f),
                ),
                f64::from_bits(0x40602cee8e01e190),
            ),
            f64::from_bits(0x4032eb254fae6dc0),
        ),
        f64::from_bits(0x3ff0000000000000),
    ))
}

#[inline]
fn inverse_norm_low<const FUSED: bool>(p: f64) -> f64 {
    pj2024_tail_from_r::<FUSED>((-p.ln()).sqrt())
}
#[inline]
fn inverse_norm_cdf_impl<const FUSED: bool>(p: f64) -> f64 {
    let u = p - 0.5;
    if u.abs() < 0.3413447460685429 {
        return inverse_norm_mid::<FUSED>(u);
    }
    if u > 0.0 {
        -inverse_norm_low::<FUSED>(1.0 - p)
    } else {
        inverse_norm_low::<FUSED>(p)
    }
}
#[inline]
pub(crate) fn inverse_log_probability(logp: f64) -> f64 {
    if logp <= -2.0 {
        pj2024_tail_from_r::<true>((-logp).sqrt())
    } else {
        inverse_norm_cdf(logp.exp())
    }
}
