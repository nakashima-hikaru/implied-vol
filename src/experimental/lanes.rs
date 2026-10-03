//! Two independent binary64 evaluation lanes with unchanged rounding points.
use core::ops::{Add, Div, Mul, Neg, Sub};

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
use core::arch::aarch64::*;

#[derive(Clone, Copy)]
pub(crate) struct F64x2(
    #[cfg(all(target_arch = "aarch64", target_feature = "neon"))] float64x2_t,
    #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))] [f64; 2],
);

impl F64x2 {
    #[inline(always)]
    pub(crate) fn new(a: f64, b: f64) -> Self {
        #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
        // This branch requires NEON at compile time; both representations have 128 bits.
        unsafe {
            Self(core::mem::transmute::<[f64; 2], float64x2_t>([a, b]))
        }
        #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
        Self([a, b])
    }
    #[inline(always)]
    pub(crate) fn splat(x: f64) -> Self {
        Self::new(x, x)
    }
    #[inline(always)]
    pub(crate) fn to_array(self) -> [f64; 2] {
        #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
        unsafe {
            core::mem::transmute::<float64x2_t, [f64; 2]>(self.0)
        }
        #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
        self.0
    }
    #[inline(always)]
    pub(crate) fn mul_add(self, b: Self, c: Self) -> Self {
        #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
        unsafe {
            Self(vfmaq_f64(c.0, self.0, b.0))
        }
        #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
        Self([
            self.0[0].mul_add(b.0[0], c.0[0]),
            self.0[1].mul_add(b.0[1], c.0[1]),
        ])
    }
}

macro_rules! binary_op {
    ($trait:ident, $method:ident, $intrinsic:ident, $op:tt) => {
        impl $trait for F64x2 {
            type Output = Self;
            #[inline(always)]
            fn $method(self, b: Self) -> Self {
                #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
                unsafe { Self($intrinsic(self.0, b.0)) }
                #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
                Self([self.0[0] $op b.0[0], self.0[1] $op b.0[1]])
            }
        }
    };
}
binary_op!(Add, add, vaddq_f64, +);
binary_op!(Sub, sub, vsubq_f64, -);
binary_op!(Mul, mul, vmulq_f64, *);
binary_op!(Div, div, vdivq_f64, /);

impl Neg for F64x2 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
        unsafe {
            Self(vnegq_f64(self.0))
        }
        #[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
        Self([-self.0[0], -self.0[1]])
    }
}
