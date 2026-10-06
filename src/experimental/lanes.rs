//! Two independent binary64 lanes, implemented entirely with `fearless_simd`.
//!
//! Use the target's guaranteed SIMD baseline; no per-operation CPU detection is
//! needed. Keep fused operations precise even on backends without hardware FMA.
use core::ops::{Add, Div, Mul, Neg, Sub};
use fearless_simd::{f64x2, prelude::*};

#[cfg(not(any(
    all(target_arch = "aarch64", target_feature = "neon"),
    all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_feature = "sse2",
        target_feature = "fxsr"
    ),
    all(target_arch = "wasm32", target_feature = "simd128")
)))]
use fearless_simd::Fallback as Backend;
#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
use fearless_simd::Neon as Backend;
#[cfg(all(
    any(target_arch = "x86", target_arch = "x86_64"),
    target_feature = "sse2",
    target_feature = "fxsr"
))]
use fearless_simd::Sse2 as Backend;
#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
use fearless_simd::WasmSimd128 as Backend;

#[inline(always)]
fn backend() -> Backend {
    #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
    return fearless_simd::Level::baseline().as_neon().unwrap();
    #[cfg(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_feature = "sse2",
        target_feature = "fxsr"
    ))]
    return fearless_simd::Level::baseline().as_sse2().unwrap();
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    return fearless_simd::Level::baseline().as_wasm_simd128().unwrap();
    #[cfg(not(any(
        all(target_arch = "aarch64", target_feature = "neon"),
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "sse2",
            target_feature = "fxsr"
        ),
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    Backend::new()
}

#[derive(Clone, Copy)]
pub(crate) struct F64x2(f64x2<Backend>);

impl F64x2 {
    #[inline(always)]
    pub(crate) fn new(a: f64, b: f64) -> Self {
        Self(f64x2::simd_from(backend(), [a, b]))
    }
    #[inline(always)]
    pub(crate) fn splat(x: f64) -> Self {
        Self(f64x2::splat(backend(), x))
    }
    #[inline(always)]
    pub(crate) fn to_array(self) -> [f64; 2] {
        self.0.into()
    }
    #[inline(always)]
    pub(crate) fn mul_add(self, b: Self, c: Self) -> Self {
        // Ordinary `mul_add` may round twice on SSE2, WASM and scalar backends.
        // Compensated Horner requires a single rounding on every platform.
        Self(self.0.mul_add_precise(b.0, c.0))
    }
}

macro_rules! binary_op {
    ($trait:ident, $method:ident, $op:tt) => {
        impl $trait for F64x2 {
            type Output = Self;
            #[inline(always)]
            fn $method(self, b: Self) -> Self {
                Self(self.0 $op b.0)
            }
        }
    };
}
binary_op!(Add, add, +);
binary_op!(Sub, sub, -);
binary_op!(Mul, mul, *);
binary_op!(Div, div, /);

impl Neg for F64x2 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::F64x2;

    fn assert_bits(actual: F64x2, expected: [f64; 2]) {
        assert_eq!(
            actual.to_array().map(f64::to_bits),
            expected.map(f64::to_bits)
        );
    }

    #[test]
    fn lanes_preserve_independent_binary64_operations() {
        let values = [
            0.0,
            -0.0,
            f64::from_bits(1),
            f64::MIN_POSITIVE,
            1.0,
            -3.0,
            f64::MAX,
            f64::INFINITY,
        ];
        for &a in &values {
            for &b in &values {
                let x = F64x2::new(a, -b);
                let y = F64x2::new(b, a);
                // NaN payloads are not part of the solver's precision contract.
                for (got, expected) in [
                    ((x + y).to_array(), [a + b, -b + a]),
                    ((x - y).to_array(), [a - b, -b - a]),
                    ((x * y).to_array(), [a * b, -b * a]),
                    ((x / y).to_array(), [a / b, -b / a]),
                ] {
                    for lane in 0..2 {
                        if expected[lane].is_nan() {
                            assert!(got[lane].is_nan());
                        } else {
                            assert_eq!(got[lane].to_bits(), expected[lane].to_bits());
                        }
                    }
                }
                assert_bits(-x, [-a, b]);
            }
            assert_bits(F64x2::splat(a), [a, a]);
        }
    }

    #[test]
    fn fused_rounding_survives_cancellation_overflow_and_underflow() {
        let e = f64::EPSILON;
        for (a, b, c) in [
            (1.0 + e, 1.0 - e, -1.0),
            (f64::MAX, 2.0, -f64::MAX),
            (f64::MIN_POSITIVE, 0.5, -f64::MIN_POSITIVE),
            (f64::from_bits(1), 0.5, f64::from_bits(1)),
            (-0.0, 1.0, -0.0),
        ] {
            let got = F64x2::new(a, c).mul_add(F64x2::new(b, b), F64x2::new(c, a));
            assert_bits(got, [a.mul_add(b, c), c.mul_add(b, a)]);
        }
        assert_eq!((1.0 + e).mul_add(1.0 - e, -1.0), -e * e);
    }
}
