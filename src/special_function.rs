mod erf_cody;
mod normal_distribution;

use crate::special_function::erf_cody::{erf_cody, erfc_cody, erfcx_cody};
use crate::special_function::normal_distribution::{erfinv, inverse_norm_cdf, norm_pdf};

pub trait SpecialFn {
    fn erf(x: f64) -> f64;
    fn erfc(x: f64) -> f64;
    fn erfcx(x: f64) -> f64;
    fn erfinv(x: f64) -> f64;
    fn inverse_norm_cdf(x: f64) -> f64;
    fn norm_pdf(x: f64) -> f64;
}

pub struct DefaultSpecialFn;
impl SpecialFn for DefaultSpecialFn {
    #[inline(always)]
    fn erf(x: f64) -> f64 {
        erf_cody(x)
    }

    #[inline(always)]
    fn erfc(x: f64) -> f64 {
        erfc_cody(x)
    }
    #[inline(always)]
    fn erfcx(x: f64) -> f64 {
        erfcx_cody(x)
    }
    #[inline(always)]
    fn erfinv(x: f64) -> f64 {
        erfinv(x)
    }
    #[inline(always)]
    fn inverse_norm_cdf(x: f64) -> f64 {
        inverse_norm_cdf(x)
    }
    #[inline(always)]
    fn norm_pdf(x: f64) -> f64 {
        norm_pdf(x)
    }
}
