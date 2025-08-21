use crate::fused_multiply_add::MulAdd;

#[inline(always)]
pub fn householder3_factor(v: f64, h2: f64, h3: f64) -> f64 {
    v.mul_add2(0.5 * h2, 1.0) / v.mul_add2(h3 / 6.0, h2).mul_add2(v, 1.0)
}

#[inline(always)]
pub fn householder4_factor(v: f64, h2: f64, h3: f64, h4: f64) -> f64 {
    v.mul_add2(h3 / 6.0, h2).mul_add2(v, 1.0)
        / v.mul_add2(h4 / 24.0, h2.mul_add2(h2 / 4.0, h3 / 3.0))
            .mul_add2(v, 1.5 * h2)
            .mul_add2(v, 1.0)
}