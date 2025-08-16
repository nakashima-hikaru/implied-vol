#[cfg(feature = "cxx_bench")]
#[cxx::bridge]
pub mod ffi {
    unsafe extern "C++" {
        include!("lets_be_rational.h");

        #[must_use] pub fn ImpliedBlackVolatility(price: f64, F: f64, K: f64, T: f64, q: f64) -> f64;
        #[must_use] pub fn Black(F: f64, K: f64, sigma: f64, T: f64, q: f64) -> f64;

        include!("ImpliedNormalVolatility.cpp");

        #[must_use] pub fn ImpliedNormalVolatility(price: f64, F: f64, K: f64, T: f64, q: f64) -> f64;
        #[must_use] pub fn Bachelier(F: f64, K: f64, sigma: f64, T: f64, q: f64) -> f64;
    }
}
