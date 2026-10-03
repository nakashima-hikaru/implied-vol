use super::*;

#[inline]
fn gap_fast(a: f64, b: f64) -> Pair {
    let r = -a / 2.0;
    let mut p = exponential_coefficients[12].hi;
    for k in (0..12).rev() {
        p = fma(p, r, exponential_coefficients[k].hi);
    }
    add(
        add(sub(Pair::from(1.0), Pair::from(b)), Pair::from(r)),
        scale(mul(Pair::from(r), Pair::from(r)), p),
    )
}
#[inline]
fn gap_refined(a: f64, b: f64) -> Pair {
    let r = -a / 2.0;
    let mut p = exponential_coefficients[18];
    for k in (0..18).rev() {
        p = add(mul(p, Pair::from(r)), exponential_coefficients[k]);
    }
    add(
        add(sub(Pair::from(1.0), Pair::from(b)), Pair::from(r)),
        mul(mul(Pair::from(r), Pair::from(r)), p),
    )
}
#[inline]
fn exp_pair(z: f64, full: bool) -> Pair {
    let k = (z * f64::from_bits(0x3ff71547652b82fe)).round_ties_even() as i32;
    let r = sub(
        Pair::from(z),
        scale(
            Pair::new(
                f64::from_bits(0x3fe62e42fefa39ef),
                f64::from_bits(0x3c7abc9e3b39803f),
            ),
            k as f64,
        ),
    );
    let mut P;
    if full {
        P = FACINV[22];
        for j in (0..22).rev() {
            P = add(mul(P, r), FACINV[j]);
        }
    } else {
        let mut p = FACINV[18].hi;
        for j in (3..18).rev() {
            p = fma(p, r.hi, FACINV[j].hi);
        }
        P = Pair::from(p);
        for j in (0..3).rev() {
            P = add(mul(P, r), FACINV[j]);
        }
    }
    Pair::new(scalbn(P.hi, k), scalbn(P.lo, k))
}
fn boundary(a: f64, b: f64) -> (f64, bool) {
    let (status, d) = crate::experimental::large::resolve_cap(a, b);
    if status != 0 {
        return (f64::NAN, status != 3 && status != 4);
    }
    let L = -d.ln();
    if !(L >= 60.0 && L <= 150.0) {
        return (f64::NAN, true);
    }
    let mut x = (2.0 * L).sqrt();
    for _ in 0..5 {
        if !(x >= 8.0 && x <= 18.0) {
            return (f64::NAN, true);
        }
        let y = fma(x, x, 2.0 * a).sqrt();
        let S = mills24_pair(x, y);
        let f = ((0.5 * x * x + f64::from_bits(0x3fed67f1c864beb5)) - S.ln()) - L;
        let derivative = (1.0 + x / y) / S;
        x -= f / derivative;
    }
    if !(x >= 8.0 && x <= 18.0) {
        return (f64::NAN, true);
    }
    (x + fma(x, x, 2.0 * a).sqrt(), true)
}
// Same continued fractions as mills24(x) + mills24(y), evaluated in parallel.
// The fixed half-open range lets LLVM expand the stages with constant numerators.
#[inline]
fn mills24_pair(x: f64, y: f64) -> f64 {
    let xy = F64x2::new(x, y);
    let mut r = F64x2::splat(0.0);
    for k in (0..24).rev() {
        r = F64x2::splat((k + 1) as f64) / (xy + r);
    }
    let q = (F64x2::splat(1.0) / (xy + r)).to_array();
    q[0] + q[1]
}
pub(super) fn small_upper(a: f64, b: f64) -> (f64, bool) {
    let mut gap = gap_fast(a, b);
    let r = -a / 2.0;
    let r2 = r * r;
    let margin = fma(
        f64::from_bits(0x3ca0000000000000),
        r2,
        f64::from_bits(0x3ba0000000000000),
    );
    if gap.hi.abs() <= margin {
        gap = gap_refined(a, b);
        if a > f64::from_bits(0x3d70000000000000)
            && gap.hi.abs() <= f64::from_bits(0x3a50000000000000)
        {
            let q = boundary(a, b);
            if q.1 {
                return q;
            }
        }
        if gap.hi.abs() <= f64::from_bits(0x39f0000000000000) * r2 {
            return (f64::NAN, false);
        }
    }
    if !(gap.hi > 0.0) {
        return (f64::NAN, false);
    }
    let A = a * a;
    let c = 1.0
        + A * (1.0 / 8.0
            + A * (1.0 / 384.0 + A * (1.0 / 46080.0 + A * (1.0 / 10321920.0 + A / 3715891200.0))));
    let p = gap.value() / (2.0 * c);
    if !(p > 0.0 && p <= f64::from_bits(0x3fc44ed0bb7cb20b)) {
        return (f64::NAN, false);
    }
    let u = -2.0 * inverse_norm_cdf(p);
    if !(u >= 2.0 && u <= 40.0) {
        return (f64::NAN, false);
    }
    let rr = A / (u * u);
    let v = u * u;
    let ratio = 1.0
        + rr * (-0.5
            + rr * ((v - 10.0) / 48.0 + rr * (-(2.0 * v * v - 23.0 * v + 246.0) / 1440.0)));
    let s = u * ratio;
    let cap = add(Pair::from(b), gap);
    let delta = div(gap, cap);
    (finish_upper(a, delta, s), false)
}
#[inline]
fn price_gap(a: f64, b: f64, cap0: f64, boundaries: bool) -> Result<(Pair, Pair), (f64, bool)> {
    let mut cap = Pair::from(cap0);
    let mut gap = sub(cap, Pair::from(b));
    if gap.hi.abs() <= f64::from_bits(0x3d20000000000000) * cap.hi {
        cap = exp_pair(-a / 2.0, false);
        gap = sub(cap, Pair::from(b));
        if gap.hi.abs() <= f64::from_bits(0x3c50000000000000) * cap.hi {
            cap = exp_pair(-a / 2.0, true);
            gap = sub(cap, Pair::from(b));
            if boundaries && gap.hi.abs() <= f64::from_bits(0x3a50000000000000) * cap.hi {
                let q = boundary(a, b);
                if q.1 {
                    return Err(q);
                }
            }
            if gap.hi.abs() <= f64::from_bits(0x39f0000000000000) * cap.hi {
                return Err((f64::NAN, false));
            }
        }
    }
    if !(gap.hi > 0.0) {
        return Err((f64::NAN, false));
    }
    Ok((cap, gap))
}
pub(super) fn wide_upper(a: f64, b: f64, cap0: f64) -> (f64, bool) {
    let (cap, gap) = match price_gap(a, b, cap0, true) {
        Ok(q) => q,
        Err(q) => return q,
    };
    let delta = div(gap, cap);
    let d = delta.value();
    if !(d > 0.0 && d <= 0.5) {
        return (f64::NAN, false);
    }
    let z = -d.ln();
    if !(z >= 0.69 && z <= 80.0) {
        return (f64::NAN, false);
    }
    let ai = if a <= 2.0 { 0 } else { 1 };
    let zi = if z < 1.4 {
        0
    } else if z < 3.0 {
        1
    } else if z < 6.0 {
        2
    } else if z < 12.0 {
        3
    } else if z < 24.0 {
        4
    } else if z < 48.0 {
        5
    } else {
        6
    };
    let k = 7 * ai + zi;
    let r = rectangles[k];
    let seed = rank_seed::<8>(factors[k], (a - r[0]) * r[1], (z - r[2]) * r[3]);
    (finish_upper(a, delta, seed), false)
}
#[inline]
fn band(delta: Pair) -> usize {
    let d = delta.hi + delta.lo.abs();
    if d <= f64::from_bits(0x3cf0000000000000) {
        5
    } else if d <= f64::from_bits(0x3df0000000000000) {
        4
    } else if d <= f64::from_bits(0x3ef0000000000000) {
        3
    } else if d <= f64::from_bits(0x3f70000000000000) {
        2
    } else if d <= f64::from_bits(0x3fb0000000000000) {
        1
    } else {
        0
    }
}
#[inline]
fn mills_sum_fast(mut x: f64, mut y: f64, band: usize) -> f64 {
    let kx = (2.0 * x).floor() as i32;
    let ky = (2.0 * y).floor() as i32;
    if !(0..32).contains(&kx) || !(0..32).contains(&ky) {
        return f64::NAN;
    }
    let a = &M_COEFF[(kx + 12) as usize];
    let b = &M_COEFF[(ky + 12) as usize];
    x -= (2 * kx + 1) as f64 * 0.25;
    y -= (2 * ky + 1) as f64 * 0.25;
    let n = upper_degrees[band][kx as usize].max(upper_degrees[band][ky as usize]);
    let h = F64x2::new(x, y);
    let mut p = F64x2::new(a[n].hi, b[n].hi);
    for j in (0..n).rev() {
        p = p.mul_add(h, F64x2::new(a[j].hi, b[j].hi));
    }
    let p = p.to_array();
    p[0] + p[1]
}
#[inline]
fn compact_mills_sum(x: f64, xl: f64, y: f64, yl: f64, band: usize) -> Pair {
    let kx = (2.0 * x).floor() as i32;
    let ky = (2.0 * y).floor() as i32;
    if !(0..32).contains(&kx) || !(0..32).contains(&ky) {
        return Pair::new(f64::NAN, f64::NAN);
    }
    let a = &M_COEFF[(kx + 12) as usize];
    let b = &M_COEFF[(ky + 12) as usize];
    let xr = x - (2 * kx + 1) as f64 * 0.25;
    let yr = y - (2 * ky + 1) as f64 * 0.25;
    let n = upper_degrees[band][kx as usize].max(upper_degrees[band][ky as usize]);
    let h = F64x2::new(xr, yr);
    let mut p = F64x2::new(a[n].hi, b[n].hi);
    for j in (1..n).rev() {
        p = p.mul_add(h, F64x2::new(a[j].hi, b[j].hi));
    }
    let prod = p * h;
    let er = p.mul_add(h, -prod);
    let c = F64x2::new(a[0].hi, b[0].hi);
    let cl = F64x2::new(a[0].lo, b[0].lo);
    let sm = prod + c;
    let v = sm - prod;
    // Keep coefficient low inside the sum residual, before product error.
    let lo = er + (((prod - (sm - v)) + (c - v)) + cl);
    let lo = F64x2::new(x, y)
        .mul_add(sm, F64x2::splat(-1.0))
        .mul_add(F64x2::new(xl, yl), lo)
        .to_array();
    let h = sm.to_array();
    let p = Pair::new(h[0], lo[0]);
    let q = Pair::new(h[1], lo[1]);
    let sum = p.hi + q.hi;
    let v = sum - p.hi;
    Pair::new(sum, ((p.hi - (sum - v)) + (q.hi - v)) + (p.lo + q.lo))
}
#[inline]
fn finish_upper(a: f64, delta: Pair, s: f64) -> f64 {
    if delta.hi + delta.lo.abs() <= f64::from_bits(0x3fb0000000000000) {
        let h = a / s;
        let t = 0.5 * s;
        let x = t - h;
        let y = t + h;
        if !(x >= 0.0 && y < 16.0 && s > 0.0 && s.is_finite()) {
            return f64::NAN;
        }
        let S = mills_sum_fast(x, y, band(delta));
        let p = f64::from_bits(0x3fd9884533d43651) * (-x * x / 2.0).exp();
        let n = (fma(p, S, -delta.hi) - delta.lo) / (s * p);
        let hh = h * h;
        let tt = t * t;
        let H2 = hh - tt;
        let H3 = fma(H2, H2, -3.0 * hh - tt);
        let ds = n * (1.0 + H2 * n / 2.0) / (1.0 + n * (H2 + H3 * n / 6.0));
        return fma(s, ds, s);
    }
    let inv = 1.0 / s;
    let h = a * inv;
    let t = 0.5 * s;
    let il = fma(-inv, s, 1.0) * inv;
    let hl = fma(a, inv, -h) + a * il;
    let x = t - h;
    let y = t + h;
    if !(x >= 0.0 && y < 16.0 && s > 0.0 && s.is_finite()) {
        return f64::NAN;
    }
    let xl = (-h - (x - t)) - hl;
    let yl = (h - (y - t)) + hl;
    let S = compact_mills_sum(x, xl, y, yl, band(delta));
    let exponent = -0.5 * x * x;
    let el = fma(-0.5 * x, x, -exponent) - x * xl;
    let ev = exponent.exp();
    let p = f64::from_bits(0x3fd9884533d43651) * ev;
    let pl = fma(f64::from_bits(0x3fd9884533d43651), ev, -p)
        + fma(p, el, -f64::from_bits(0x3c7cbc0d30ebfd15) * ev);
    let mut residual = fma(p, S.hi, -delta.hi);
    residual = fma(p, S.lo, residual);
    residual = fma(pl, S.hi, residual);
    residual -= delta.lo;
    let n = residual / (s * p);
    let hh = h * h;
    let tt = t * t;
    let H2 = hh - tt;
    let H3 = fma(H2, H2, -3.0 * hh - tt);
    let ds = n * (1.0 + H2 * n / 2.0) / (1.0 + n * (H2 + H3 * n / 6.0));
    fma(s, ds, s)
}
pub(super) fn shared(a: f64, b: f64) -> f64 {
    let (cap, gap) = match price_gap(a, b, (-a / 2.0).exp(), false) {
        Ok(q) => q,
        Err(_) => return f64::NAN,
    };
    let delta = div(gap, cap);
    let d = delta.value();
    if !(d > 0.0 && d <= 0.5) {
        return f64::NAN;
    }
    let mut x = -inverse_norm_cdf(d / 2.0);
    let v = x * x;
    let w = 1.0 / (1.0 + v);
    let r = x / (v + 2.0 * a).sqrt();
    if !(w <= 0.69 && r >= 0.13) {
        return f64::NAN;
    }
    let xx = fma(w, 2.0 / 0.69, -1.0);
    let yy = (r - 0.565) / 0.435;
    let correction = tensor_seed::<11, 11>(&GAP_PACKED, xx, yy);
    x = (v + correction).sqrt();
    let y = (x * x + 2.0 * a).sqrt();
    let S = dual_mills(Pair::from(x), Pair::from(y), true).value();
    let p = f64::from_bits(0x3fd9884533d43651) * (-x * x / 2.0).exp();
    let n = ((fma(p, S, -delta.hi) - delta.lo) / p) / (1.0 + x / y);
    let l = -x + (y - x) / (y * y);
    let lp = -1.0 + (x / y - 1.0) / (y * y) - 2.0 * x * (y - x) / (y * y * y * y);
    let H3 = fma(l, l, lp);
    let dx = n * (1.0 + l * n / 2.0) / (1.0 + n * (l + H3 * n / 6.0));
    let X = add(Pair::from(x), Pair::from(dx));
    let Y = add(mul(X, X), Pair::from(2.0 * a));
    let yy = Y.hi.sqrt();
    let y = Pair::new(yy, (fma(-yy, yy, Y.hi) + Y.lo) / (2.0 * yy));
    add(X, y).value()
}
