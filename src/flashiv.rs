//! `FlashIV` log-price Householder inversion.
//!
//! Independently implemented from Le Floc'h and Healy, arXiv:2605.29102v1,
//! equations (logc), (logvega), (d2), (d3), (h3), and Algorithm 1:
//! <https://arxiv.org/abs/2605.29102v1>
//! Li's seed coefficients are from Li and Lee, MPRA 6867, equation (42):
//! <https://mpra.ub.uni-muenchen.de/6867/1/MPRA_paper_6867.pdf>
//! Abramowitz--Stegun 7.1.26 supplies the inexpensive erfcx polynomial.
//!
//! The pure entry follows Algorithm 1: a terminal microscopic-price guard;
//! otherwise one cheap H3 step and two exact H3 steps, with a third exact step
//! when the residual entering the second exact step is at least 1e-4. The
//! upper-price guard uses three Halley steps on the complementary log price.
//! The available May 2026 author sources use two upper Householder steps;
//! Algorithm 1 takes precedence here. Their shared microscopic guard supplies
//! the Mills-ratio and local Black--Bachelier formulas and the LFK2026 rational
//! seed coefficients (the `FlashIV` text cites the older 2016 approximation):
//! <https://chasethedevil.github.io/post/thiophene-iv-rust-full.zip>
//!
//! Prices use sqrt-forward beta space: the x/2 terms in the log-price objective
//! cancel. Special functions and opt-in FMA follow this crate's provider and
//! arithmetic policy, so this is not a bit-identical port of the author code.
//! There is no adaptive convergence loop, inverse-method fallback, or final
//! `FlashIV+` price alignment. Arithmetic guards retain the current iterate,
//! as in the available author implementation; the microscopic terminal guard
//! may return its defensive zero-volatility limit.
//!
//! The default hybrid's restricted entry retains its qualified scaled Black
//! expansions in the exact steps. The pure entry uses the paper's direct
//! erfcx log objective instead, including its cancellation limits.

// Preserve the crate's opt-in FMA policy.
#![allow(clippy::suboptimal_flops)]

use crate::fused_multiply_add::MulAdd;
use crate::lets_be_rational::bs_option_price::scaled_normalised_black_and_ln_vega;
use crate::lets_be_rational::bs_option_price::uses_scaled_expansion;
use crate::lets_be_rational::special_function::SpecialFn;
use std::f64::consts::FRAC_1_SQRT_2;
use std::f64::consts::LN_2;

const SQRT_2_OVER_PI: f64 = 0.797_884_560_802_865_4;
const INV_SQRT_PI: f64 = 0.564_189_583_547_756_3;
const LN_2_PI: f64 = 1.837_877_066_409_345_3;
#[cfg(feature = "flashiv")]
const SQRT_2_PI: f64 = 2.506_628_274_631_000_2;
#[cfg(feature = "flashiv")]
const UPPER_PRICE_GUARD: f64 = f64::from_bits(0.99_f64.to_bits() - 1);

/// Pure Algorithm 1 inversion, without the optional final price correction.
#[cfg(feature = "flashiv")]
#[inline]
pub fn normalised_volatility<SpFn: SpecialFn>(beta: f64, theta_x: f64, b_max: f64) -> Option<f64> {
    if !(theta_x <= 0.0 && theta_x.is_finite() && beta > 0.0 && beta < b_max && b_max.is_finite()) {
        return None;
    }
    let m = -theta_x;
    let c = beta / b_max;
    let ln_beta = beta.ln();
    if c <= 1e-6 && m <= 1e-8 {
        return Some(microscopic_volatility::<SpFn>(m, beta, ln_beta));
    }
    let ln_c = ln_beta + 0.5 * m;
    let mut v = paper_initial_guess(theta_x, c, ln_c).max(1e-10);
    if c >= UPPER_PRICE_GUARD {
        v = v.max(upper_price_seed(theta_x, c));
        let ln_gap = (b_max - beta).ln();
        v = upper_halley_step::<SpFn>(theta_x, v, ln_gap);
        v = upper_halley_step::<SpFn>(theta_x, v, ln_gap);
        v = upper_halley_step::<SpFn>(theta_x, v, ln_gap);
    } else {
        v = paper_h3_step::<SpFn, true>(theta_x, v, ln_beta).0;
        v = paper_h3_step::<SpFn, false>(theta_x, v, ln_beta).0;
        let (next, residual_on_second_entry) = paper_h3_step::<SpFn, false>(theta_x, v, ln_beta);
        v = next;
        if residual_on_second_entry.abs() >= 1e-4 {
            v = paper_h3_step::<SpFn, false>(theta_x, v, ln_beta).0;
        }
    }
    (v > 0.0 && v.is_finite()).then_some(v)
}

#[cfg(feature = "flashiv")]
#[inline(always)]
fn paper_initial_guess(x: f64, c: f64, ln_c: f64) -> f64 {
    let m = -x;
    if m < 3.0 && c > 0.0005 && c < UPPER_PRICE_GUARD {
        li_seed(x, c)
    } else if c <= 0.0005 && m < 0.01 {
        x.mul_add2(x, (SQRT_2_PI * c) * (SQRT_2_PI * c)).sqrt()
    } else if c >= UPPER_PRICE_GUARD {
        upper_price_seed(x, c)
    } else if c <= 0.5 && ln_c < -2.0 {
        let d2 = -2.0 * ln_c - LN_2_PI;
        (-2.0 * x) / (d2.sqrt() + (d2 - 2.0 * x).sqrt())
    } else {
        // The put-side seed follows the available author implementation.
        let put = (-x).exp() - c;
        if put > 1e-300 {
            let d2 = -2.0 * put.ln() - 2.0 * x - LN_2_PI;
            if d2 > 0.0 {
                return d2.sqrt() + (d2 + 2.0 * x).max(0.0).sqrt();
            }
        }
        (2.0 * m).sqrt()
    }
}

#[cfg(feature = "flashiv")]
#[inline(always)]
fn upper_price_seed(x: f64, c: f64) -> f64 {
    let r = -2.0 * (-c).ln_1p() - LN_2_PI;
    let root = r.sqrt();
    let d = root - (0.5 * r.ln()) * root / (r + 1.0);
    d + (d * d - 2.0 * x).max(0.0).sqrt()
}

// Pure FlashIV always evaluates the direct erfcx difference. Keep this helper
// separate from the default hybrid's qualified scaled-expansion objective.
#[cfg(feature = "flashiv")]
#[inline(always)]
fn paper_h3_step<SpFn: SpecialFn, const FAST: bool>(x: f64, v: f64, ln_beta: f64) -> (f64, f64) {
    let h = x / v;
    let t = 0.5 * v;
    let delta = if FAST {
        fast_erfcx(-FRAC_1_SQRT_2 * (h + t)) - fast_erfcx(-FRAC_1_SQRT_2 * (h - t))
    } else {
        SpFn::erfcx(-FRAC_1_SQRT_2 * (h + t)) - SpFn::erfcx(-FRAC_1_SQRT_2 * (h - t))
    };
    if !(delta > 0.0 && delta.is_finite()) {
        return (v, f64::INFINITY);
    }
    let h2 = h * h;
    let t2 = t * t;
    let log_beta = if !FAST && x == 0.0 {
        SpFn::erf(FRAC_1_SQRT_2 * t).ln()
    } else {
        -0.5 * (h2 + t2) - LN_2 + delta.ln()
    };
    let residual = log_beta - ln_beta;
    let d1 = SQRT_2_OVER_PI / delta;
    let a = h2 - t2;
    let vov = (h + t) * (h - t) / v;
    let d2 = vov - d1;
    let d3 = (a.mul_add2(a, (-3.0).mul_add2(h2, -t2))) / (v * v) - 3.0 * d1 * vov + 2.0 * d1 * d1;
    let nu = -residual / d1;
    let numerator = nu.mul_add2(0.5 * d2, 1.0);
    let denominator = nu.mul_add2(nu * (d3 / 6.0), nu.mul_add2(d2, 1.0));
    let next = v + nu * (numerator / denominator);
    (
        if denominator != 0.0 && next > 0.0 && next.is_finite() {
            next
        } else {
            v
        },
        residual,
    )
}

#[cfg(feature = "flashiv")]
#[inline(always)]
fn upper_halley_step<SpFn: SpecialFn>(x: f64, v: f64, ln_beta_gap: f64) -> f64 {
    let h = x / v;
    let d1 = h + 0.5 * v;
    let z1 = FRAC_1_SQRT_2 * d1;
    let z2 = FRAC_1_SQRT_2 * (v - d1);
    let sum = SpFn::erfcx(z1) + SpFn::erfcx(z2);
    if !(sum > 0.0 && sum.is_finite()) {
        return v;
    }
    // ln(beta_max-beta(v)) = -d1^2/2 + x/2 - ln(2) + ln(sum).
    let residual = -z1 * z1 + 0.5 * x - LN_2 + sum.ln() - ln_beta_gap;
    let inverse_slope = (0.5 * SQRT_2_PI) * sum;
    let newton = inverse_slope * residual;
    let d1_prime = 0.5 - h / v;
    let correction = residual * (d1 * d1_prime * inverse_slope - 1.0);
    let denominator = 1.0 - 0.5 * correction;
    // The published shared guard retains a Newton step when the Halley
    // denominator is too small, and retains v for an invalid update.
    let step = if denominator > 0.1 {
        newton / denominator
    } else {
        newton
    };
    let next = v + step;
    if next > 0.0 && next.is_finite() {
        next
    } else {
        v
    }
}

#[cfg(feature = "flashiv")]
fn microscopic_volatility<SpFn: SpecialFn>(m: f64, beta: f64, ln_beta: f64) -> f64 {
    if let Some(v) = bachelier_mills_tail::<SpFn>(m, ln_beta) {
        return v;
    }
    let mut v = bachelier_rational_seed(m, beta);
    if !(v > 0.0 && v.is_finite()) {
        return 0.0;
    }
    if m > 0.0 && m / v > 4.0 {
        return v;
    }
    for _ in 0..2 {
        let a = m / v;
        let phi = (-0.5 * a * a).exp() / SQRT_2_PI;
        let derivative = phi * (-0.125 * v * v).exp();
        if !(derivative > 0.0 && derivative.is_finite()) {
            break;
        }
        let i0 = v * phi - m * (0.5 * SpFn::erfc(FRAC_1_SQRT_2 * a));
        let v2 = v * v;
        let m2 = m * m;
        let i2 = (v2 * v * phi - m2 * i0) / 3.0;
        let i4 = (v2 * v2 * v * phi - m2 * i2) / 5.0;
        let estimate = i0 - 0.125 * i2 + i4 / 128.0;
        let step = (estimate - beta) / derivative;
        let next = v - step;
        if !(next > 0.0 && next.is_finite()) {
            break;
        }
        v = next;
        if step.abs() <= 2.0 * (v.next_up() - v) {
            break;
        }
    }
    v
}

#[cfg(feature = "flashiv")]
fn bachelier_mills_tail<SpFn: SpecialFn>(m: f64, ln_beta: f64) -> Option<f64> {
    if !(m > 0.0 && ln_beta.is_finite()) {
        return None;
    }
    let log_ratio = ln_beta - m.ln();
    if log_ratio >= -20.0 {
        return None;
    }
    let mut a = (-2.0 * (log_ratio + 0.5 * LN_2_PI)).max(1.0).sqrt();
    for _ in 0..4 {
        let mills = (0.5 * std::f64::consts::PI).sqrt() * SpFn::erfcx(FRAC_1_SQRT_2 * a);
        let gap = a.recip() - mills;
        if !(gap > 0.0 && gap.is_finite()) {
            return None;
        }
        let residual = -0.5 * a * a - 0.5 * LN_2_PI + gap.ln() - log_ratio;
        let slope = -a + (1.0 - (a * a).recip() - a * mills) / gap;
        let step = residual / slope;
        let next = a - step;
        if !(next > 0.0 && next.is_finite()) {
            return None;
        }
        a = next;
        if step.abs() <= 1e-14 * a.max(1.0) {
            break;
        }
    }
    let v = m / a;
    (v.is_finite() && m / v > 4.0).then_some(v)
}

// Numerical coefficients in the author's May 2026 shared guard (LFK2026),
// not the historical LFK2016 approximation cited by the FlashIV paper. The
// piecewise rational basis is v/(m+beta) in the central region and v/m in
// three transformed log-price regions. Coefficients are highest degree first.
#[cfg(feature = "flashiv")]
fn bachelier_rational_seed(m: f64, beta: f64) -> f64 {
    if m == 0.0 {
        return SQRT_2_PI * beta;
    }
    let z = beta / m;
    if z > 0.15 {
        return (beta + m) * rational(m / beta, &BACHELIER_CENTRAL_P, &BACHELIER_CENTRAL_Q);
    }
    if z < 1e-300 {
        return 0.0;
    }
    let eta =
        -(z.ln() + 1.897_119_984_885_881_3) / (690.775_527_898_213_7 - 1.897_119_984_885_881_3);
    let (p, q) = if eta < 0.01 {
        (&BACHELIER_TAIL1_P, &BACHELIER_TAIL1_Q)
    } else if eta < 0.09 {
        (&BACHELIER_TAIL2_P, &BACHELIER_TAIL2_Q)
    } else {
        (&BACHELIER_TAIL3_P, &BACHELIER_TAIL3_Q)
    };
    m * rational(eta, p, q)
}

#[cfg(feature = "flashiv")]
#[inline(always)]
fn rational(x: f64, numerator: &[f64], denominator: &[f64]) -> f64 {
    let horner = |coefficients: &[f64]| {
        coefficients[1..]
            .iter()
            .fold(coefficients[0], |v, &c| v.mul_add2(x, c))
    };
    horner(numerator) / horner(denominator)
}

#[cfg(feature = "flashiv")]
const BACHELIER_CENTRAL_P: [f64; 10] = [
    -2.707_970_425_512_152e-09,
    5.038_038_374_003_797e-06,
    0.000_587_928_493_313_025_6,
    0.018_395_449_137_629_65,
    0.236_012_991_833_862_85,
    1.471_719_297_921_368_2,
    4.815_324_741_041_999,
    8.387_836_492_959_392,
    7.320_714_511_547_501,
    2.506_628_274_631_000_2,
];

#[cfg(feature = "flashiv")]
const BACHELIER_CENTRAL_Q: [f64; 9] = [
    1.067_153_385_375_751_3e-05,
    0.000_908_750_187_647_611_8,
    0.022_318_863_699_157_96,
    0.230_574_991_341_749_2,
    1.175_343_674_707_846_7,
    3.181_652_961_948_590_3,
    4.636_111_360_386_343_5,
    3.420_542_541_404_568_5,
    1.0,
];

#[cfg(feature = "flashiv")]
const BACHELIER_TAIL1_P: [f64; 12] = [
    -3.483_088_522_889_665_5e+18,
    1.102_993_700_042_91e+18,
    -4.024_137_161_683_780_5e+17,
    -7.437_143_499_458_291e+16,
    -2_494_778_633_228_567.5,
    -28_533_126_098_770.652,
    -51_081_746_906.337_14,
    2_405_330_731.308_365_3,
    30_786_109.258_591_59,
    195_717.866_058_198_88,
    682.178_299_087_715_9,
    1.249_105_594_466_411_1,
];

#[cfg(feature = "flashiv")]
const BACHELIER_TAIL1_Q: [f64; 10] = [
    -8.658_682_189_272_106e+18,
    -5.436_913_992_018_88e+17,
    -9_163_742_300_293_810.0,
    -47_328_368_624_205.15,
    446_201_521_087.863_83,
    9_146_846_117.682_764,
    71_809_151.130_195_66,
    322_084.448_691_854_1,
    831.825_318_775_392_4,
    1.0,
];

#[cfg(feature = "flashiv")]
const BACHELIER_TAIL2_P: [f64; 12] = [
    2_991_471_363_897.162_6,
    -8_766_255_362_542.491,
    30_942_671_490_695.523,
    57_382_549_836_842.99,
    20_380_699_276_606.562,
    2_580_761_776_362.431,
    139_249_460_078.530_36,
    3_450_744_505.068_248_3,
    41_266_059.989_788_875,
    257_450.335_132_154_6,
    882.143_141_269_435_3,
    1.248_604_954_973_179_2,
];

#[cfg(feature = "flashiv")]
const BACHELIER_TAIL2_Q: [f64; 10] = [
    2_107_197_868_545_195.2,
    1_455_887_995_322_584.8,
    292_823_593_402_229.5,
    23_120_099_325_362.99,
    795_484_765_720.863_6,
    12_479_445_436.712_564,
    95_610_288.258_178_26,
    417_679.926_778_023_36,
    991.194_200_989_211_3,
    1.0,
];

#[cfg(feature = "flashiv")]
const BACHELIER_TAIL3_P: [f64; 12] = [
    -101.865_471_163_877_16,
    2_995.245_780_148_021,
    -104_999.322_668_089_05,
    -1_904_380.098_369_818_2,
    -6_493_743.087_133_431,
    -7_655_631.945_038_427,
    -3_620_913.058_956_987_7,
    -678_815.730_118_479_5,
    -36_372.743_382_534_16,
    1_536.691_127_704_749_3,
    134.977_888_075_358_55,
    1.619_476_154_121_739_7,
];

#[cfg(feature = "flashiv")]
const BACHELIER_TAIL3_Q: [f64; 10] = [
    -22_419_688.389_442_544,
    -151_588_857.602_825_4,
    -293_510_878.971_269_7,
    -216_291_778.562_989_83,
    -64_430_925.073_288_23,
    -6_854_155.825_363_629_5,
    -34_483.659_597_472_05,
    22_752.425_724_960_64,
    591.576_136_257_250_4,
    1.0,
];

#[inline]
pub fn try_lowest_branch<SpFn: SpecialFn>(beta: f64, theta_x: f64, b_max: f64) -> Option<f64> {
    if !(theta_x <= -0.01 && theta_x.is_finite() && beta > 0.0 && b_max > 0.0) {
        return None;
    }
    let c = beta / b_max;
    if !(c > 0.0 && c < 0.5) {
        return None;
    }
    let ln_beta = beta.ln();
    let ln_c = ln_beta - 0.5 * theta_x;
    let seed = initial_guess(theta_x, c, ln_c).max(1e-10);
    let (v1, _) = h3_step::<SpFn, true>(theta_x, seed, ln_beta)?;
    let (v2, _) = h3_step::<SpFn, false>(theta_x, v1, ln_beta)?;
    let (v3, residual_on_second_entry) = h3_step::<SpFn, false>(theta_x, v2, ln_beta)?;
    if residual_on_second_entry.abs() >= 1e-4 {
        h3_step::<SpFn, false>(theta_x, v3, ln_beta).map(|(v, _)| v)
    } else {
        Some(v3)
    }
}

#[inline(always)]
fn initial_guess(x: f64, c: f64, ln_c: f64) -> f64 {
    let abs_x = -x;
    if abs_x < 3.0 && c > 0.0005 {
        li_seed(x, c)
    } else if c <= 0.5 && ln_c < -2.0 {
        let d2 = -2.0 * ln_c - LN_2_PI;
        (-2.0 * x) / (d2.sqrt() + (d2 - 2.0 * x).sqrt())
    } else {
        (2.0 * abs_x).sqrt()
    }
}

#[inline(always)]
fn li_seed(x: f64, c: f64) -> f64 {
    // Collect sum_{i+j<=3} a_ij x^i c^j as Horner polynomials in c,
    // with coefficient polynomials in x, respecting the opt-in FMA choice.
    let m0 = x
        .mul_add2(-0.216_197_632_156_68, 0.089_753_944_048_51)
        .mul_add2(x, -0.406_619_903_654_27)
        .mul_add2(x, -0.000_061_030_981_65);
    let m1 = x
        .mul_add2(3.838_158_853_945_65, -36.194_052_215_990_28)
        .mul_add2(x, 5.339_676_433_576_88);
    let m2 = x.mul_add2(41.217_726_327_328_34, 3.250_234_253_323_60);
    let numerator = c
        .mul_add2(83.845_932_244_177_96, m2)
        .mul_add2(c, m1)
        .mul_add2(c, m0);
    let n0 = x
        .mul_add2(-0.033_269_442_900_44, 0.430_276_195_531_68)
        .mul_add2(x, -0.484_665_363_616_20)
        .mul_add2(x, 1.0);
    let n1 = x
        .mul_add2(-0.047_638_023_588_53, -1.341_022_799_820_50)
        .mul_add2(x, 22.963_021_090_107_94);
    let n2 = x.mul_add2(2.457_825_742_942_44, -0.772_688_245_324_68);
    let denominator = c
        .mul_add2(-5.705_315_006_451_09, n2)
        .mul_add2(c, n1)
        .mul_add2(c, n0);
    numerator / denominator
}

#[inline(always)]
fn h3_step<SpFn: SpecialFn, const FAST: bool>(
    x: f64,
    v: f64,
    ln_beta_target: f64,
) -> Option<(f64, f64)> {
    if !(v > 0.0 && v.is_finite()) {
        return None;
    }
    let inv_v = v.recip();
    let h = x * inv_v;
    let t = 0.5 * v;
    let (log_beta, d1) = if FAST || !uses_scaled_expansion(h, t) {
        let z_plus = -FRAC_1_SQRT_2 * (h + t);
        let z_minus = -FRAC_1_SQRT_2 * (h - t);
        let delta = if FAST {
            fast_erfcx(z_plus) - fast_erfcx(z_minus)
        } else {
            SpFn::erfcx(z_plus) - SpFn::erfcx(z_minus)
        };
        if !(delta > 0.0 && delta.is_finite()) {
            return None;
        }
        (
            (-0.5).mul_add2(t.mul_add2(t, h * h), delta.ln() - LN_2),
            SQRT_2_OVER_PI / delta,
        )
    } else {
        let (scaled_b, ln_vega) = scaled_normalised_black_and_ln_vega::<SpFn>(0.5 * x, h, t);
        if !(scaled_b > 0.0 && scaled_b.is_finite() && ln_vega.is_finite()) {
            return None;
        }
        (scaled_b.ln() + ln_vega, scaled_b.recip())
    };
    let residual = log_beta - ln_beta_target;
    let h2 = h * h;
    let t2 = t * t;
    let h2_minus_t2 = h2 - t2;
    let b2 = h2_minus_t2 * inv_v;
    let d2 = b2 - d1;
    let b3 = h2_minus_t2.mul_add2(h2_minus_t2, (-3.0).mul_add2(h2, -t2)) * (inv_v * inv_v);
    let d3 = (-3.0 * d1).mul_add2(d2, b3 - d1 * d1);
    let nu = -residual / d1;
    let numerator = nu.mul_add2(0.5 * d2, 1.0);
    let denominator = nu.mul_add2(nu.mul_add2(d3 / 6.0, d2), 1.0);
    let next = (nu * (numerator / denominator)) + v;
    if next > 0.0 && next.is_finite() {
        Some((next, residual))
    } else {
        None
    }
}

#[inline(always)]
fn fast_erfcx(z: f64) -> f64 {
    if z < 0.0 {
        return (z * z).exp().mul_add2(2.0, -fast_erfcx(-z));
    }
    if z >= 2.5 {
        let inv_z = z.recip();
        let w = inv_z * inv_z;
        let correction = w.mul_add2(-1.875, 0.75).mul_add2(w, -0.5).mul_add2(w, 1.0);
        return (INV_SQRT_PI * inv_z) * correction;
    }
    let t = z.mul_add2(0.327_591_1, 1.0).recip();
    t.mul_add2(1.061_405_429, -1.453_152_027)
        .mul_add2(t, 1.421_413_741)
        .mul_add2(t, -0.284_496_736)
        .mul_add2(t, 0.254_829_592)
        * t
}
