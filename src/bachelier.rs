use crate::constants::{ONE_OVER_SQRT_TWO_PI, SQRT_TWO_PI};
use crate::normal_distribution::norm_pdf;
use std::cmp::Ordering;

#[inline(always)]
fn intrinsic_value(forward: f64, strike: f64, q: bool) -> f64 {
    (if !q {
        strike - forward
    } else {
        forward - strike
    })
    .max(0.0)
    .abs()
}

#[inline(always)]
fn phi_tilde_times_x(x: f64) -> f64 {
    if x.abs() <= 0.612_003_180_962_480_7 {
        let h = (x * x - 1.872_739_467_540_974_8E-1) * 5.339_771_053_755_08;
        let g = (1.964_154_984_377_470_3E-1
            + h * (2.944_481_222_626_891_4E-3 + 3.095_828_855_856_471E-5 * h))
            / (1.0
                + h * (3.026_101_684_659_232_6E-2
                    + h * (3.373_546_191_189_62E-4
                        + h * (1.290_112_376_540_573_2E-6 - 1.671_197_583_524_420_5E-9 * h))));
        return ONE_OVER_SQRT_TWO_PI + x * (0.5 + x * g);
    }

    if x.is_sign_positive() {
        return phi_tilde_times_x(-x) + x;
    }

    if x >= -3.5 {
        let g = (3.989_422_804_009_617_5E-1
            + x * ((-2.882_725_012_271_64E-1)
                + x * (1.174_893_477_005_507_4E-1
                    + x * ((-2.920_893_049_832_423_4E-2)
                        + x * (4.670_481_708_734_893E-3
                            + x * ((-4.444_840_548_247_636E-4)
                                + x * (1.986_526_744_238_593_6E-5
                                    + x * (7.638_739_347_414_361E-10
                                        + 1.329_152_522_013_758_3E-11 * x))))))))
            / (1.0
                + x * ((-1.975_906_139_672_860_6)
                    + x * (1.770_933_219_893_362_5
                        + x * ((-9.435_025_002_644_624E-1)
                            + x * (3.281_611_814_538_859_5E-1
                                + x * ((-7.669_740_808_821_474E-2)
                                    + x * (1.184_322_430_309_622_3E-2
                                        + x * ((-1.115_141_636_552_486_1E-3)
                                            + 4.974_100_533_375_869E-5 * x))))))));
        return (-0.5 * (x * x)).exp() * g;
    }

    let w = (x * x).recip();
    let g = (2.999_999_999_999_991
        + w * (2.365_455_662_782_315E2
            + w * (6.812_677_344_935_879E3
                + w * (8.969_794_159_836_079E4
                    + w * (5.516_392_059_126_862E5
                        + w * (1.434_506_112_333_566_2E6
                            + w * (1.150_498_824_634_488_2E6 + 1.186_760_040_099_769_1E4 * w)))))))
        / (1.0
            + w * (8.384_852_209_273_714E1
                + w * (2.655_135_058_780_958E3
                    + w * (4.055_529_088_467_379E4
                        + w * (3.166_737_476_299_376_6E5
                            + w * (1.232_979_595_802_432_2E6
                                + w * (2.140_981_054_061_905E6 + 1.214_566_780_409_316E6 * w)))))));
    ONE_OVER_SQRT_TWO_PI * (-0.5 * (x * x)).exp() * w * (1.0 - g * w)
}

#[inline(always)]
fn phi_tilde(x: f64) -> f64 {
    phi_tilde_times_x(x) / x
}

#[inline(always)]
fn inv_phi_tilde(phi_tilde_star: f64) -> f64 {
    if phi_tilde_star > 1.0 {
        return -inv_phi_tilde(1.0 - phi_tilde_star);
    }
    if !phi_tilde_star.is_sign_negative() {
        return f64::NAN;
    }
    let x_bar = if phi_tilde_star < -0.00188203927 {
        // Equation (2.1)
        let g = (phi_tilde_star - 0.5).recip();
        let g2 = g * g;
        // Equation (2.2)
        let xi_bar = (0.032114372355
            - g2 * (0.016969777977 - g2 * (0.002620733246 - 0.000096066952861 * g2)))
            / (1.0 - g2 * (0.6635646938 - g2 * (0.14528712196 - 0.010472855461 * g2)));
        // Equation (2.3)
        g * (ONE_OVER_SQRT_TWO_PI + xi_bar * g2)
    } else {
        // Equation (2.4)
        let h = (-(-phi_tilde_star).ln()).sqrt();
        // Equation (2.5)
        (9.4883409779 - h * (9.6320903635 - h * (0.58556997323 + 2.1464093351 * h)))
            / (1.0 - h * (0.65174820867 + h * (1.5120247828 + 0.000066437847132 * h)))
    };
    // Equation (2.7)
    let q = (phi_tilde(x_bar) - phi_tilde_star) / norm_pdf(x_bar);
    let x2 = x_bar * x_bar;
    // Equation (2.6)
    x_bar
        + 3.0 * q * x2 * (2.0 - q * x_bar * (2.0 + x2))
            / (6.0
                + q * x_bar * (-12.0 + x_bar * (6.0 * q + x_bar * (-6.0 + q * x_bar * (3.0 + x2)))))
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
pub(crate) fn bachelier(forward: f64, strike: f64, sigma: f64, t: f64, q: bool) -> f64 {
    let s = sigma.abs() * t.sqrt();
    if s < f64::MIN_POSITIVE {
        return intrinsic_value(forward, strike, q);
    }
    let moneyness = if q {
        forward - strike
    } else {
        strike - forward
    };
    let x = moneyness / s;
    s * phi_tilde_times_x(x)
}

#[inline(always)]
pub(crate) fn implied_normal_volatility(
    price: f64,
    forward: f64,
    strike: f64,
    t: f64,
    q: bool,
) -> f64 {
    if forward == strike {
        return price * SQRT_TWO_PI / t.sqrt();
    }
    let intrinsic = intrinsic_value(forward, strike, q);
    match price.total_cmp(&intrinsic) {
        Ordering::Less => f64::NEG_INFINITY,
        Ordering::Equal => 0.0,
        Ordering::Greater => {
            let absolute_moneyness = (forward - strike).abs();
            let phi_tilde_star = (intrinsic - price) / absolute_moneyness;
            let x_star = inv_phi_tilde(phi_tilde_star);
            absolute_moneyness / (x_star * t.sqrt()).abs()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn reconstruction_call_atm() {
        for i in 1..100 {
            let price = 0.01 * i as f64;
            let f = 100.0;
            let k = f;
            let t = 1.0;
            let q = true;
            let sigma = implied_normal_volatility(price, f, k, t, q);
            let reprice = bachelier(f, k, sigma, t, q);
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
            let q = false;
            let sigma = implied_normal_volatility(price, f, k, t, q);
            let reprice = bachelier(f, k, sigma, t, q);
            assert!((price - reprice).abs() < 5e-14);
        }
    }

    #[test]
    fn reconstruction_random_call_intrinsic() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.gen();
            let price = 1e5 * r2;
            let f = r + 1e5 * r2;
            let k = f - price;
            let t = 1e5 * r3;
            let q = true;
            let sigma = implied_normal_volatility(price, f, k, t, q);
            let reprice = bachelier(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_call_itm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.gen();
            let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
            let f = 1.0;
            let k = 1.0 * r;
            let t = 1e5 * r3;
            let q = true;
            let sigma = implied_normal_volatility(price, f, k, t, q);
            let reprice = bachelier(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_call_otm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.gen();
            let price = 1.0 * r * r2;
            let f = 1.0 * r;
            let k = 1.0;
            let t = 1e5 * r3;
            let q = true;
            let sigma = implied_normal_volatility(price, f, k, t, q);
            let reprice = bachelier(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_put_itm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.gen();
            let price = 1.0 * r * r2;
            let f = 1.0;
            let k = 1.0 * r;
            let t = 1e5 * r3;
            let q = false;
            let sigma = implied_normal_volatility(price, f, k, t, q);
            let reprice = bachelier(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }

    #[test]
    fn reconstruction_random_put_otm() {
        let n = 100_000;
        let seed: [u8; 32] = [13; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        for _ in 0..n {
            let (r, r2, r3): (f64, f64, f64) = rng.gen();
            let price = 1.0 * (1.0 - r) + 1.0 * r * r2;
            let f = 1.0 * r;
            let k = 1.0;
            let t = 1e5 * r3;
            let q = false;
            let sigma = implied_normal_volatility(price, f, k, t, q);
            let reprice = bachelier(f, k, sigma, t, q);
            assert!((price - reprice).abs() <= 2.0 * f64::EPSILON);
        }
    }
}
