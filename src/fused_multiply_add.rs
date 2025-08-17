pub trait MulAdd {
    fn mul_add2(self, a: f64, b: f64) -> f64;
}

impl MulAdd for f64 {
    #[cfg(not(feature = "fma"))]
    #[inline(always)]
    fn mul_add2(self, a: f64, b: f64) -> f64 {
        a * self + b
    }

    #[cfg(feature = "fma")]
    #[inline(always)]
    fn mul_add2(self, a: f64, b: f64) -> f64 {
        self.mul_add(a, b)
    }
}
