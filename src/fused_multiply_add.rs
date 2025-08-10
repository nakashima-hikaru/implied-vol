pub(super) trait MulAdd {
    fn mul_add2(self, a: f64, b: f64) -> f64;
}

#[cfg(feature = "fma")]
#[inline]
fn fma_function(r: f64, a: f64, b: f64) -> f64 {
    r.mul_add(a, b)
}

#[cfg(not(feature = "fma"))]
#[inline]
fn fma_function(r: f64, a: f64, b: f64) -> f64 {
    a * r + b
}

impl MulAdd for f64 {
    #[inline]
    fn mul_add2(self, a: f64, b: f64) -> f64 {
        fma_function(self, a, b)
    }
}
