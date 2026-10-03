//! v40 dispatcher retained by experimental for a <= 10, with independent SIMD evaluation lanes.
#![allow(non_snake_case)]
mod central;
mod constants;
mod upper;
use crate::experimental::lanes::F64x2;
use crate::experimental::math::inverse_norm_cdf;
use constants::*;

#[derive(Clone, Copy, Debug)]
struct Pair {
    hi: f64,
    lo: f64,
}
impl Pair {
    const fn new(hi: f64, lo: f64) -> Self {
        Self { hi, lo }
    }
    const fn from(hi: f64) -> Self {
        Self { hi, lo: 0.0 }
    }
    #[inline]
    fn value(self) -> f64 {
        self.hi + self.lo
    }
}
#[inline]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    a.mul_add(b, c)
}
#[inline]
fn add(a: Pair, b: Pair) -> Pair {
    let s = a.hi + b.hi;
    let z = s - a.hi;
    let e = ((a.hi - (s - z)) + (b.hi - z)) + (a.lo + b.lo);
    let h = s + e;
    Pair::new(h, (s - h) + e)
}
#[inline]
fn sub(a: Pair, b: Pair) -> Pair {
    add(a, Pair::new(-b.hi, -b.lo))
}
#[inline]
fn mul(a: Pair, b: Pair) -> Pair {
    let h = a.hi * b.hi;
    let l = fma(a.hi, b.hi, -h);
    let l = fma(a.hi, b.lo, l);
    let l = fma(a.lo, b.hi, l);
    add(Pair::from(h), Pair::from(l + a.lo * b.lo))
}
#[inline]
fn div(a: Pair, b: Pair) -> Pair {
    let h = a.hi / b.hi;
    let r = sub(a, mul(Pair::from(h), b));
    add(Pair::from(h), Pair::from((r.hi + r.lo) / b.hi))
}
#[inline]
fn scale(a: Pair, b: f64) -> Pair {
    mul(a, Pair::from(b))
}

#[inline]
fn mills_index(x: f64) -> Option<usize> {
    let k = (2.0 * x).floor() as i32;
    if (-12..32).contains(&k) {
        Some((k + 12) as usize)
    } else {
        None
    }
}
#[inline]
fn single_mills(xx: Pair, n: usize) -> Pair {
    let Some(k) = mills_index(xx.hi) else {
        return Pair::new(f64::NAN, f64::NAN);
    };
    let c = &M_COEFF[k];
    let x = sub(xx, Pair::from((2.0 * (k as f64 - 12.0) + 1.0) * 0.25));
    let mut p = c[n].hi;
    for i in (3..n).rev() {
        p = fma(p, x.hi, c[i].hi);
    }
    let mut lo = 0.0;
    for i in (0..3).rev() {
        let prod = p * x.hi;
        let mut err = fma(p, x.hi, -prod) + p * x.lo;
        let sm = prod + c[i].hi;
        let z = sm - prod;
        err += ((prod - (sm - z)) + (c[i].hi - z)) + c[i].lo;
        lo = fma(lo, x.hi, err);
        p = sm;
    }
    Pair::new(p, lo)
}
#[inline]
fn dual_mills(x: Pair, y: Pair, sum: bool) -> Pair {
    let (Some(kx), Some(ky)) = (mills_index(x.hi), mills_index(y.hi)) else {
        return Pair::new(f64::NAN, f64::NAN);
    };
    let n = mills_degree[kx].max(mills_degree[ky]);
    let a = &M_COEFF[kx];
    let b = &M_COEFF[ky];
    let x = sub(x, Pair::from((2.0 * (kx as f64 - 12.0) + 1.0) * 0.25));
    let y = sub(y, Pair::from((2.0 * (ky as f64 - 12.0) + 1.0) * 0.25));
    let h = F64x2::new(x.hi, y.hi);
    let l = F64x2::new(x.lo, y.lo);
    let mut p = F64x2::new(a[n].hi, b[n].hi);
    for i in (3..n).rev() {
        p = p.mul_add(h, F64x2::new(a[i].hi, b[i].hi));
    }
    let mut lo = F64x2::splat(0.0);
    for i in (0..3).rev() {
        let c = F64x2::new(a[i].hi, b[i].hi);
        let cl = F64x2::new(a[i].lo, b[i].lo);
        let prod = p * h;
        let err = p.mul_add(h, -prod) + p * l;
        let sm = prod + c;
        let z = sm - prod;
        let err = err + (((prod - (sm - z)) + (c - z)) + cl);
        lo = lo.mul_add(h, err);
        p = sm;
    }
    let p = p.to_array();
    let lo = lo.to_array();
    let p0 = Pair::new(p[0], lo[0]);
    let p1 = Pair::new(p[1], lo[1]);
    if sum { add(p0, p1) } else { sub(p0, p1) }
}
#[inline]
fn dpoly(h: Pair, t: Pair) -> Pair {
    let Some(k) = mills_index(h.hi) else {
        return Pair::new(f64::NAN, f64::NAN);
    };
    let c = &M_COEFF[k];
    let x = sub(h, Pair::from((2.0 * (k as f64 - 12.0) + 1.0) * 0.25));
    let tt = mul(t, t);
    let mut e = c[18].hi;
    let mut o = 0.0;
    for i in (3..18).rev() {
        let ne = fma(e, x.hi, fma(o, tt.hi, c[i].hi));
        o = fma(o, x.hi, e);
        e = ne;
    }
    let mut E = Pair::from(e);
    let mut O = Pair::from(o);
    for i in (0..3).rev() {
        let ne = add(add(mul(E, x), mul(O, tt)), c[i]);
        O = add(E, mul(O, x));
        E = ne;
    }
    scale(O, -1.0)
}
#[inline]
fn forward_D(h: Pair, t: Pair) -> Pair {
    if t.hi <= 0.001 {
        let Some(k) = mills_index(h.hi) else {
            return Pair::new(f64::NAN, f64::NAN);
        };
        let d = sub(Pair::from(1.0), mul(h, single_mills(h, mills_degree[k])));
        let H = h.hi * h.hi;
        let T = t.hi * t.hi;
        let i3 = fma(H + 3.0, d.hi, -1.0);
        let i5 = fma(fma(H, H + 10.0, 15.0), d.hi, -H - 7.0);
        return add(d, Pair::from(T * fma(T, i5 / 120.0, i3 / 6.0)));
    }
    if t.hi <= 0.0625 {
        return dpoly(h, t);
    }
    div(dual_mills(sub(h, t), add(h, t), false), scale(t, 2.0))
}
#[inline]
fn finish(a: f64, b: f64, s: f64) -> f64 {
    let hv = a / s;
    let h = Pair::new(hv, fma(-hv, s, a) / s);
    let t = Pair::from(s * 0.5);
    if !(h.hi >= 0.0 && h.hi <= 8.0 && t.hi > 0.0 && t.hi < 5.5 && h.hi + t.hi < 15.5) {
        return f64::NAN;
    }
    let hh = mul(h, h);
    let tt = mul(t, t);
    let D = forward_D(h, t);
    if !(D.hi > 0.0 && D.hi.is_finite()) {
        return f64::NAN;
    }
    let exponent = scale(add(hh, tt), -0.5);
    let ev = exponent.hi.exp();
    let v = mul(
        mul(Pair::from(ev), Pair::new(1.0, exponent.lo)),
        Pair::new(
            f64::from_bits(0x3fd9884533d43651),
            -f64::from_bits(0x3c7cbc0d30ebfd15),
        ),
    );
    let sv = scale(v, s);
    let residual = (fma(-sv.hi, D.hi, b) - sv.hi * D.lo) - sv.lo * D.hi;
    let n = residual / sv.hi;
    let H2 = hh.hi - tt.hi;
    let H3 = fma(H2, H2, -3.0 * hh.hi - tt.hi);
    let ds = n * (1.0 + 0.5 * H2 * n) / (1.0 + n * (H2 + H3 * n / 6.0));
    fma(s, ds, s)
}
#[inline]
fn tensor_seed<const NJ: usize, const NI: usize>(c: &[[f64; 12]; NJ], x: f64, y: f64) -> f64 {
    let mut rows = [0.0; 12];
    for i in 0..NI {
        rows[i] = c[NJ - 1][i];
        for j in (0..NJ - 1).rev() {
            rows[i] = fma(rows[i], y, c[j][i]);
        }
    }
    let mut r = 0.0;
    for i in (0..NI).rev() {
        r = fma(r, x, rows[i]);
    }
    r
}
#[inline]
fn rank_seed<const DEGREE: usize>(t: (&[f64], usize), x: f64, y: f64) -> f64 {
    // The tables have four fixed shapes. Check the selected slice once so
    // Horner stages can use constant offsets without per-coefficient checks.
    const { assert!(DEGREE == 8 || DEGREE == 10) };
    let (c, rank) = t;
    if rank <= 2 {
        if DEGREE == 10 {
            rank_seed_shape::<10, 4, 44>(
                c.first_chunk::<44>().expect("rank coefficients"),
                rank,
                x,
                y,
            )
        } else {
            rank_seed_shape::<8, 4, 36>(
                c.first_chunk::<36>().expect("rank coefficients"),
                rank,
                x,
                y,
            )
        }
    } else if DEGREE == 10 {
        rank_seed_shape::<10, 8, 88>(
            c.first_chunk::<88>().expect("rank coefficients"),
            rank,
            x,
            y,
        )
    } else {
        rank_seed_shape::<8, 8, 72>(
            c.first_chunk::<72>().expect("rank coefficients"),
            rank,
            x,
            y,
        )
    }
}
#[inline]
fn rank_seed_shape<const DEGREE: usize, const LANES: usize, const LENGTH: usize>(
    c: &[f64; LENGTH],
    rank: usize,
    x: f64,
    y: f64,
) -> f64 {
    let half = LANES / 2;
    assert!(rank <= half, "rank factors");
    let mut values = [0.0; 8];
    for k in (0..LANES).step_by(2) {
        let z = F64x2::splat(if k < half { x } else { y });
        // Preserve the original degree-eight y factors in wider degree-ten seeds.
        let n = if DEGREE == 10 && LANES == 8 && k >= half {
            8
        } else {
            DEGREE
        };
        let mut v = F64x2::new(c[n * LANES + k], c[n * LANES + k + 1]);
        for j in (0..n).rev() {
            v = v.mul_add(z, F64x2::new(c[j * LANES + k], c[j * LANES + k + 1]));
        }
        values[k..k + 2].copy_from_slice(&v.to_array());
    }
    let mut v = values[0] * values[half];
    for k in 1..rank {
        v = fma(values[k], values[half + k], v);
    }
    v
}
#[inline]
fn wing_raw(a: f64, b: f64) -> f64 {
    if !(a >= 1e-100 && a <= 10.0 && b > 0.0 && b.is_finite()) {
        return f64::NAN;
    }
    let L = 2.0 * ((a / b).ln() - f64::from_bits(0x3fed67f1c864beb5));
    if !(L >= 16.0 && L <= 80.0) {
        return f64::NAN;
    }
    let k = if L < 32.0 { 0 } else { 1 };
    let x = (L - if k == 1 { 56.0 } else { 24.0 }) / if k == 1 { 24.0 } else { 8.0 };
    let y = fma(a, a / 50.0, -1.0);
    let factor = tensor_seed::<7, 9>(&WING_PACKED[k], x, y);
    finish(a, b, a / (L * factor).sqrt())
}
#[inline]
fn small_rank(a: f64, b: f64) -> f64 {
    let q = f64::from_bits(0x40040d931ff62706) * fma(0.5, a, b);
    if !(q >= 1e-100 && q <= 2.0) {
        return f64::NAN;
    }
    let w = ((a * 0.5) / b).ln_1p();
    if !(w >= 0.0 && w <= 8.5) {
        return f64::NAN;
    }
    let k = (if w < 2.0 {
        0
    } else if w < 5.0 {
        1
    } else {
        2
    }) * 2
        + if q <= 1.0 { 0 } else { 1 };
    let r = small_rects[k];
    finish(
        a,
        b,
        q * rank_seed::<10>(small_lr[k], (w - r[0]) * r[1], fma(q, q, -r[2]) * r[3]),
    )
}
#[inline]
fn finite_rank(a: f64, b: f64) -> f64 {
    let z = -b.ln() - a / 2.0;
    if !(z >= 0.35 && z <= 13.0) {
        return f64::NAN;
    }
    let k = (if a < 0.25 {
        0
    } else if a < 1.0 {
        1
    } else if a < 4.0 {
        2
    } else {
        3
    }) * 5
        + (if z < 1.2 {
            0
        } else if z < 3.0 {
            1
        } else if z < 6.0 {
            2
        } else if z < 10.0 {
            3
        } else {
            4
        });
    let r = core_rects[k];
    finish(
        a,
        b,
        rank_seed::<10>(core_lr[k], (a - r[0]) * r[1], (z - r[2]) * r[3]),
    )
}
#[inline]
fn near_atm(a: f64, b: f64) -> f64 {
    let p = add(Pair::from(b), Pair::from(a / 2.0));
    let e = crate::experimental::math::erfinv(p.hi);
    let s = mul(
        Pair::from(e),
        Pair::new(
            f64::from_bits(0x4006a09e667f3bcd),
            -f64::from_bits(0x3cabdd3413b26456),
        ),
    );
    let v = f64::from_bits(0x3fd9884533d43651) * (-s.hi * s.hi / 8.0).exp();
    let w = scale(mul(s, s), 0.125);
    let mut tail = ATM_SERIES[18].hi;
    for n in (3..18).rev() {
        tail = fma(tail, w.hi, ATM_SERIES[n].hi);
    }
    let mut P = Pair::from(tail);
    for n in (0..3).rev() {
        P = add(mul(P, w), ATM_SERIES[n]);
    }
    let price = mul(
        mul(s, P),
        Pair::new(
            f64::from_bits(0x3fd9884533d43651),
            -f64::from_bits(0x3c7cbc0d30ebfd15),
        ),
    );
    let correction = sub(p, price).value() / v;
    add(s, Pair::from(correction)).value()
}

pub(crate) fn solve(a: f64, b: f64) -> f64 {
    if !(a >= 0.0 && a <= 10.0 && a.is_finite() && b > 0.0 && b <= 1.0 && b.is_finite()) {
        return f64::NAN;
    }
    if a == 0.0 {
        if b < 1e-100 {
            return scaled_small(a, b);
        }
        let s = if b > f64::from_bits(0x3e30000000000000) {
            central::atm(b)
        } else {
            f64::from_bits(0x4006a09e667f3bcd) * crate::experimental::math::erfinv(b)
        };
        return if s > 0.0 && s < 1e300 { s } else { f64::NAN };
    }
    if a >= 1e-100 && b < a * f64::from_bits(0x3cb2200d83be3a2f) {
        let s = central::balanced(a, b);
        if positive(s) {
            return s;
        }
    }
    if a <= 0.36 && b < 0.21 {
        let s = central::deferred(a, b);
        if positive(s) {
            return s;
        }
    }
    if a <= 0.5 && b > 0.4 {
        let (s, terminal) = upper::small_upper(a, b);
        if terminal || positive(s) {
            return s;
        }
    }
    if a <= f64::from_bits(0x3e10000000000000) * b && b >= 1e-100 && b <= 0.75 {
        let s = near_atm(a, b);
        if positive(s) {
            return s;
        }
    }
    if b < a * 0.000134 {
        let s = wing_raw(a, b);
        if positive(s) {
            return s;
        }
    }
    if a >= 0.5 {
        let cap = (-a / 2.0).exp();
        if b > 0.5 * cap {
            let (s, terminal) = upper::wide_upper(a, b, cap);
            if terminal || positive(s) {
                return s;
            }
            let s = upper::shared(a, b);
            if positive(s) {
                return s;
            }
        }
    }
    if a <= 1.6 && fma(0.5, a, b) <= f64::from_bits(0x3fe9884533d43651) {
        let s = small_rank(a, b);
        if positive(s) {
            return s;
        }
    }
    if a >= 0.0625 {
        let s = finite_rank(a, b);
        if positive(s) {
            return s;
        }
    }
    if a.max(b) < 1e-100 {
        return scaled_small(a, b);
    }
    // Preserve the archived final safety fallback, implemented entirely in Rust.
    let s = crate::experimental::lbr_fallback::solve(a, b);
    if s > 0.0 && s < 1e300 && s.is_finite() {
        s
    } else {
        f64::NAN
    }
}
#[inline]
fn positive(s: f64) -> bool {
    s > 0.0 && s.is_finite()
}
fn scaled_small(a: f64, b: f64) -> f64 {
    let e = ilogb(a.max(b));
    let shift = -40 - e;
    let result = solve(scalbn(a, shift), scalbn(b, shift));
    let s = scalbn(result, -shift);
    if positive(s) { s } else { f64::NAN }
}
#[inline]
fn ilogb(x: f64) -> i32 {
    let bits = x.to_bits();
    let exponent = ((bits >> 52) & 2047) as i32;
    if exponent == 0 {
        63 - bits.leading_zeros() as i32 - 1074
    } else {
        exponent - 1023
    }
}
#[inline]
fn scalbn(mut x: f64, mut e: i32) -> f64 {
    while e > 1023 {
        x *= f64::from_bits(2046u64 << 52);
        e -= 1023;
    }
    // Keep the intermediate normal, so a subnormal result is rounded only once.
    while e < -1022 {
        x *= f64::from_bits(54u64 << 52);
        e += 969;
    }
    x * f64::from_bits(((e + 1023) as u64) << 52)
}

#[cfg(test)]
mod tests {
    use super::solve;
    #[test]
    fn archived_cpp_branch_fixtures() {
        // Expected bit patterns come from the unchanged C++ experimental with its AVX2 lane order.
        let cases: &[(u64, u64, u64)] = &[
            (0x0000000000000000, 0x3fb999999999999a, 0x3fd015abc78e92d1),
            (0x0000000000000000, 0x3feccccccccccccd, 0x400a515209676abe),
            (0x0000000000000000, 0x3cd0000000000000, 0x3ce40d931ff62706),
            (0x0000000000000000, 0x0000000000000383, 0x00000000000008cd),
            (0x16687e92154ef7ac, 0x124212dd4de70913, 0x16364f16a8097a31),
            (0x3fb999999999999a, 0x3fa999999999999a, 0x3fcd6496956a4ac3),
            (0x3fb999999999999a, 0x39b4484bfeebc2a0, 0x3f82e778a810e85c),
            (0x3ff0000000000000, 0x3ddb7cdfd9d7bdbb, 0x3fc61cc36e0fd08a),
            (0x3fb999999999999a, 0x3fd3333333333333, 0x3fecd3902c5566af),
            (0x4000000000000000, 0x3f847ae147ae147b, 0x3ff0b6994bd2a0ab),
            (0x3fc999999999999a, 0x3fe6666666666666, 0x400444dcaef09557),
            (0x4000000000000000, 0x3fd6666666666666, 0x4012bb1c928d54ff),
            (0x3fe0000000000000, 0x3fe8ebef9eac820a, 0x40309ba6b580a27b),
        ];
        for &(a, b, want) in cases {
            assert_eq!(
                solve(f64::from_bits(a), f64::from_bits(b)).to_bits(),
                want,
                "a={a:x}, b={b:x}"
            );
        }
    }
    #[test]
    fn invalid_domain() {
        for (a, b) in [
            (-1.0, 0.1),
            (0.0, 0.0),
            (0.0, 1.0),
            (1.0, 1.0),
            (f64::NAN, 0.1),
            (1.0, f64::INFINITY),
        ] {
            assert!(solve(a, b).is_nan());
        }
    }

    #[test]
    fn fallback_keeps_attainable_accuracy_for_tiny_moneyness() {
        // Direct high-precision Black root from the FMA experiment. Sharing
        // the fused seed with raw LBR raised rho from 0.306 to 10.001 here.
        let a = f64::from_bits(0x2d4e4f223fe361d6);
        let b = f64::from_bits(0x2ca4bee704377548);
        let hi = f64::from_bits(0x2d37c1936defffbd);
        let lo = f64::from_bits(0xa9dfbacb310da52e);
        let eta = f64::from_bits(0x3cb1c62650f174c3);
        let s = crate::experimental::lbr_fallback::solve(a, b);
        let rho = ((s - hi) - lo).abs() / (hi * eta);
        assert!(rho < 1.0, "fallback rho={rho}");
    }
}
