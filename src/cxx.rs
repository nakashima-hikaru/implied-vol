#[cfg(feature = "bench")]
#[cxx::bridge]
pub mod ffi {
    unsafe extern "C++" {
        include!("lets_be_rational.h");

        pub fn ImpliedBlackVolatility(price: f64, F: f64, K: f64, T: f64, q: f64) -> f64;
    }
}
