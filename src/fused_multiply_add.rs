pub trait MulAdd {
    fn mul_add2(self, a: f64, b: f64) -> f64;
}

pub trait MulAddSelected {
    fn mul_add2_sel(self, a: f64, b: f64) -> f64;
}

#[cfg(not(feature = "fma"))]
#[inline]
const fn mul_add_fast(x: f64, a: f64, b: f64) -> f64 {
    x * a + b
}

impl MulAdd for f64 {
    #[cfg(not(feature = "fma"))]
    #[inline]
    fn mul_add2(self, a: f64, b: f64) -> f64 {
        mul_add_fast(self, a, b)
    }

    #[cfg(feature = "fma")]
    #[inline]
    fn mul_add2(self, a: f64, b: f64) -> f64 {
        self.mul_add(a, b)
    }
}

impl MulAddSelected for f64 {
    #[inline]
    fn mul_add2_sel(self, a: f64, b: f64) -> f64 {
        self.mul_add2(a, b)
    }
}

#[inline]
pub fn mul_add_selected(x: f64, a: f64, b: f64) -> f64 {
    x.mul_add2(a, b)
}
