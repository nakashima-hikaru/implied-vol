use crate::fused_multiply_add::MulAdd;
use crate::lets_be_rational::constants::{FRAC_1_SQRT_2_PI, SQRT_2_PI};
use crate::lets_be_rational::special_function::SpecialFn;
use crate::lets_be_rational::special_function::normal_distribution::inv_norm_pdf;
use std::cmp::Ordering;

#[inline(always)]
fn intrinsic_value<const IS_CALL: bool>(forward: f64, strike: f64) -> f64 {
    if IS_CALL {
        forward - strike
    } else {
        strike - forward
    }
    .max(0.0)
}

#[inline(always)]
fn phi_tilde_times_x(x: f64) -> f64 {
    if x.abs() <= 0.612_003_180_962_480_7 {
        let h = x.mul_add2(x, -1.872_739_467_540_974_8E-1) * 5.339_771_053_755_08;
        let g = h
            .mul_add2(3.095_828_855_856_471E-5, 2.944_481_222_626_891_4E-3)
            .mul_add2(h, 1.964_154_984_377_470_3E-1)
            / h.mul_add2(-1.671_197_583_524_420_5E-9, 1.290_112_376_540_573_2E-6)
                .mul_add2(h, 3.373_546_191_189_62E-4)
                .mul_add2(h, 3.026_101_684_659_232_6E-2)
                .mul_add2(h, 1.0);
        return x.mul_add2(g, 0.5).mul_add2(x, FRAC_1_SQRT_2_PI);
    }

    if x > 0.0 {
        return phi_tilde_times_x(-x) + x;
    }

    if x >= -3.5 {
        let g = x
            .mul_add2(1.329_152_522_013_758_3E-11, 7.638_739_347_414_361E-10)
            .mul_add2(x, 1.986_526_744_238_593_6E-5)
            .mul_add2(x, -4.444_840_548_247_636E-4)
            .mul_add2(x, 4.670_481_708_734_893E-3)
            .mul_add2(x, -2.920_893_049_832_423_4E-2)
            .mul_add2(x, 1.174_893_477_005_507_4E-1)
            .mul_add2(x, -2.882_725_012_271_64E-1)
            .mul_add2(x, 3.989_422_804_009_617_5E-1)
            / x.mul_add2(4.974_100_533_375_869E-5, -1.115_141_636_552_486_1E-3)
                .mul_add2(x, 1.184_322_430_309_622_3E-2)
                .mul_add2(x, -7.669_740_808_821_474E-2)
                .mul_add2(x, 3.281_611_814_538_859_5E-1)
                .mul_add2(x, -9.435_025_002_644_624E-1)
                .mul_add2(x, 1.770_933_219_893_362_5)
                .mul_add2(x, -1.975_906_139_672_860_6)
                .mul_add2(x, 1.0);
        return (-0.5 * x * x).exp() * g;
    }

    let w = (x * x).recip();
    let g = w
        .mul_add2(1.186_760_040_099_769_1E4, 1.150_498_824_634_488_2E6)
        .mul_add2(w, 1.434_506_112_333_566_2E6)
        .mul_add2(w, 5.516_392_059_126_862E5)
        .mul_add2(w, 8.969_794_159_836_079E4)
        .mul_add2(w, 6.812_677_344_935_879E3)
        .mul_add2(w, 2.365_455_662_782_315E2)
        .mul_add2(w, 2.999_999_999_999_991)
        / w.mul_add2(1.214_566_780_409_316E6, 2.140_981_054_061_905E6)
            .mul_add2(w, 1.232_979_595_802_432_2E6)
            .mul_add2(w, 3.166_737_476_299_376_6E5)
            .mul_add2(w, 4.055_529_088_467_379E4)
            .mul_add2(w, 2.655_135_058_780_958E3)
            .mul_add2(w, 8.384_852_209_273_714E1)
            .mul_add2(w, 1.0);

    FRAC_1_SQRT_2_PI * (-0.5 * x * x).exp() * w * g.mul_add2(-w, 1.0)
}

#[inline(always)]
fn phi_tilde(x: f64) -> f64 {
    phi_tilde_times_x(x) / x
}

#[inline(always)]
fn inv_phi_tilde<SpFn: SpecialFn>(phi_tilde_star: f64) -> f64 {
    if phi_tilde_star > 1.0 {
        return -inv_phi_tilde::<SpFn>(1.0 - phi_tilde_star);
    }
    let x_bar = if phi_tilde_star < -0.001_882_039_27 {
        // Equation (2.1)
        let g = (phi_tilde_star - 0.5).recip();
        let g2 = g * g;
        // Equation (2.2)
        let xi_bar = g2.mul_add2(
            -g2.mul_add2(
                -0.000_096_066_952_861_f64.mul_add2(-g2, 0.002_620_733_246),
                0.016_969_777_977,
            ),
            0.032_114_372_355,
        ) / g2
            .mul_add2(-0.010_472_855_461, 0.145_287_121_96)
            .mul_add2(-g2, 0.663_564_693_8)
            .mul_add2(-g2, 1.0);
        // Equation (2.3)
        g * xi_bar.mul_add2(g2, FRAC_1_SQRT_2_PI)
    } else {
        // Equation (2.4)
        let h = (-(-phi_tilde_star).ln()).sqrt();
        // Equation (2.5)
        h.mul_add2(2.146_409_335_1, 0.585_569_973_23)
            .mul_add2(-h, 9.632_090_363_5)
            .mul_add2(-h, 9.488_340_977_9)
            / h.mul_add2(0.000_066_437_847_132, 1.512_024_782_8)
                .mul_add2(h, 0.651_748_208_67)
                .mul_add2(-h, 1.0)
    };
    // Equation (2.7)
    let q = (phi_tilde(x_bar) - phi_tilde_star) * inv_norm_pdf(x_bar);
    let x2 = x_bar * x_bar;
    // Equation (2.6)
    x_bar
        + 3.0 * q * x2 * (q * x_bar).mul_add2(-(2.0 + x2), 2.0)
            / (q * x_bar).mul_add2(
                x_bar.mul_add2(
                    6.0f64.mul_add2(q, x_bar * (q * x_bar).mul_add2(3.0 + x2, -6.0)),
                    -12.0,
                ),
                6.0,
            )
}

/// Calculates the price of an option using Bachelier's model.
///
/// # Arguments
///
/// * `forward` - The forward price of the underlying asset.
/// * `strike` - The strike price of the option.
/// * `sigma` - The volatility of the underlying asset.
/// * `t` - The time to expiration of the option.
/// * `q` - A boolean flag indicating whether the option is a put (true) or a call (false).
///
/// # Returns
///
/// The price of the option.
pub fn bachelier_price<const IS_CALL: bool>(forward: f64, strike: f64, sigma: f64, t: f64) -> f64 {
    assert!(!forward.is_nan() && !strike.is_nan() && sigma >= 0.0 && t >= 0.0);
    let s = sigma * t.sqrt();
    if s == 0.0 {
        return intrinsic_value::<IS_CALL>(forward, strike);
    }
    let moneyness = if IS_CALL {
        forward - strike
    } else {
        strike - forward
    };
    let x = moneyness / s;
    s * phi_tilde_times_x(x)
}

#[inline(always)]
pub fn implied_normal_volatility_input_unchecked<SpFn: SpecialFn, const IS_CALL: bool>(
    price: f64,
    forward: f64,
    strike: f64,
    t: f64,
) -> Option<f64> {
    if forward == strike {
        return Some(price * SQRT_2_PI / t.sqrt());
    }
    let intrinsic = intrinsic_value::<IS_CALL>(forward, strike);
    match price.total_cmp(&intrinsic) {
        Ordering::Less => None,
        Ordering::Equal => Some(0.0),
        Ordering::Greater => {
            let absolute_moneyness = (forward - strike).abs();
            let phi_tilde_star = (intrinsic - price) / absolute_moneyness;
            let x_star = inv_phi_tilde::<SpFn>(phi_tilde_star);
            Some(absolute_moneyness / (x_star * t.sqrt()).abs())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lets_be_rational::special_function::DefaultSpecialFn;
    use rand::Rng;

    #[test]
    fn reconstruction_call_atm() {
        for i in 1..100 {
            let price = 0.01 * i as f64;
            let f = 100.0;
            let k = f;
            let t = 1.0;
            let sigma =
                implied_normal_volatility_input_unchecked::<DefaultSpecialFn, true>(price, f, k, t)
                    .unwrap();
            let reprice = bachelier_price::<true>(f, k, sigma, t);
            assert!((price - reprice).abs() < 5e-14);
        }
    }

    #[test]
    fn reconstruction_put_atm() {
        for i in 1..100 {
            let price = 0.01 * i as f64;
            let f = 100.0;
            let k = f;
            let t = 1.0;
            let sigma = implied_normal_volatility_input_unchecked::<DefaultSpecialFn, false>(
                price, f, k, t,
            )
            .unwrap();
            let reprice = bachelier_price::<false>(f, k, sigma, t);
            assert!((price - reprice).abs() < 5e-14);
        }
    }

    #[test]
    fn reconstruction_random_call_intrinsic() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1e5 * r2;
            let f = r + 1e5 * r2;
            let k = f - price;
            let t = 1e5 * r3;
            let sigma =
                implied_normal_volatility_input_unchecked::<DefaultSpecialFn, true>(price, f, k, t)
                    .unwrap();
            let reprice = bachelier_price::<true>(f, k, sigma, t);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_call_itm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
            let f = 1.0;
            let k = 1.0 * r;
            let t = 1e5 * r3;
            let sigma =
                implied_normal_volatility_input_unchecked::<DefaultSpecialFn, true>(price, f, k, t)
                    .unwrap();
            let reprice = bachelier_price::<true>(f, k, sigma, t);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_call_otm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * r * r2;
            let f = 1.0 * r;
            let k = 1.0;
            let t = 1e5 * r3;
            let sigma =
                implied_normal_volatility_input_unchecked::<DefaultSpecialFn, true>(price, f, k, t)
                    .unwrap();
            let reprice = bachelier_price::<true>(f, k, sigma, t);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_put_itm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * r * r2;
            let f = 1.0;
            let k = 1.0 * r;
            let t = 1e5 * r3;
            let sigma = implied_normal_volatility_input_unchecked::<DefaultSpecialFn, false>(
                price, f, k, t,
            )
            .unwrap();
            let reprice = bachelier_price::<false>(f, k, sigma, t);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_put_otm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.random();
            let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
            let f = 1.0 * r;
            let k = 1.0;
            let t = 1e5 * r3;
            let sigma = implied_normal_volatility_input_unchecked::<DefaultSpecialFn, false>(
                price, f, k, t,
            )
            .unwrap();
            let reprice = bachelier_price::<false>(f, k, sigma, t);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }
}
