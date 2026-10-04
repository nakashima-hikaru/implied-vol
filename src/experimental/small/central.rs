use super::*;

#[inline]
fn rawsum(a: Pair, b: Pair) -> Pair {
    let s = a.hi + b.hi;
    let v = s - a.hi;
    Pair::new(s, ((a.hi - (s - v)) + (b.hi - v)) + (a.lo + b.lo))
}
#[inline]
fn rawmul(a: Pair, b: Pair) -> Pair {
    let p = a.hi * b.hi;
    Pair::new(p, fma(a.hi, b.hi, -p) + a.hi * b.lo + a.lo * b.hi)
}
#[inline]
fn rawdiv(a: Pair, b: Pair) -> Pair {
    let q = a.hi / b.hi;
    Pair::new(q, ((fma(-q, b.hi, a.hi) + a.lo) - q * b.lo) / b.hi)
}
#[inline]
fn renorm(p: Pair) -> Pair {
    let s = p.hi + p.lo;
    Pair::new(s, p.lo - (s - p.hi))
}
#[inline]
fn step(a: Pair, x: Pair, c: Pair) -> Pair {
    let p = a.hi * x.hi;
    let ep = fma(a.hi, x.hi, -p);
    let s = p + c.hi;
    let z = s - p;
    let es = (p - (s - z)) + (c.hi - z);
    let e = fma(a.lo, x.hi, fma(a.hi, x.lo, (ep + es) + c.lo));
    Pair::new(s, e)
}
#[inline]
fn sinhc_half(a: f64) -> Pair {
    let x = 0.5 * a;
    let wh = x * x;
    let w = Pair::new(wh, fma(x, x, -wh));
    let mut t = 1.0 / 355687428096000.0;
    for c in [
        1.0 / 1307674368000.0,
        1.0 / 6227020800.0,
        1.0 / 39916800.0,
        1.0 / 362880.0,
    ] {
        t = fma(t, wh, c);
    }
    let mut p = Pair::from(t);
    for c in [
        Pair::new(
            f64::from_bits(0x3f2a01a01a01a01a),
            f64::from_bits(0x3b6a01a01a01a01a),
        ),
        Pair::new(
            f64::from_bits(0x3f81111111111111),
            f64::from_bits(0x3c01111111111111),
        ),
        Pair::new(
            f64::from_bits(0x3fc5555555555555),
            f64::from_bits(0x3c65555555555555),
        ),
        Pair::from(1.0),
    ] {
        p = step(p, w, c);
    }
    renorm(p)
}
#[inline]
fn normal_deferred(z: Pair) -> Pair {
    let k = if z.hi < 0.16 {
        0
    } else if z.hi < 0.30 {
        1
    } else if z.hi < 0.40 {
        2
    } else {
        3
    };
    const MID: [Pair; 4] = [
        Pair::new(
            f64::from_bits(0x3fb47ae147ae147b),
            -f64::from_bits(0x3c3eb851eb851eb8),
        ),
        Pair::new(
            f64::from_bits(0x3fcd70a3d70a3d71),
            -f64::from_bits(0x3c670a3d70a3d70a),
        ),
        Pair::new(
            f64::from_bits(0x3fd6666666666666),
            f64::from_bits(0x3c7999999999999a),
        ),
        Pair::new(
            f64::from_bits(0x3fdccccccccccccd),
            -f64::from_bits(0x3c6999999999999a),
        ),
    ];
    const SCALE: [Pair; 4] = [
        Pair::from(12.5),
        Pair::new(
            f64::from_bits(0x402c924924924925),
            -f64::from_bits(0x3ccb6db6db6db6db),
        ),
        Pair::from(20.0),
        Pair::from(20.0),
    ];
    let t = rawmul(rawsum(z, Pair::new(-MID[k].hi, -MID[k].lo)), SCALE[k]);
    // Independent Horner chains share the argument; keep each chain's degree
    // and every original rounding point while evaluating in two NEON lanes.
    let mut value_derivative = F64x2::new(fma(0.0, t.hi, NC[k][25]), 25.0 * NC[k][25]);
    let arg = F64x2::splat(t.hi);
    for i in (2..25).rev() {
        value_derivative = value_derivative.mul_add(arg, F64x2::new(NC[k][i], i as f64 * NC[k][i]));
    }
    let value_derivative = value_derivative.to_array();
    let dp = fma(value_derivative[1], t.hi, NC[k][1]);
    let mut p = Pair::from(value_derivative[0]);
    for i in (0..2).rev() {
        let prod = p.hi * t.hi;
        let e = fma(p.hi, t.hi, -prod);
        let s = prod + NC[k][i];
        let residual = prod - (s - NC[k][i]);
        p = Pair::new(s, fma(p.lo, t.hi, (e + residual) + NL[k][i]));
    }
    Pair::new(p.hi, fma(dp, t.lo, p.lo))
}
#[inline]
fn conversion_pair<const LEFT: usize, const RIGHT: usize>(
    a: &[f64; LEFT],
    b: &[f64; RIGHT],
    x: f64,
) -> [f64; 2] {
    // The rows have different degrees. Evaluate the longer row's leading
    // stages first, then pair the remaining stages without padding either row.
    const { assert!(LEFT >= RIGHT && RIGHT > 0) };
    let head = LEFT - RIGHT;
    let mut p = a[0];
    for c in &a[1..=head] {
        p = fma(p, x, *c);
    }
    let mut p = F64x2::new(p, b[0]);
    let arg = F64x2::splat(x);
    for i in 1..RIGHT {
        p = p.mul_add(arg, F64x2::new(a[head + i], b[i]));
    }
    p.to_array()
}

pub(super) fn deferred(a: f64, b: f64) -> f64 {
    // These margins necessarily fail the q or z range guard below:
    // b<a/16 gives z>0.503; b+a/2>=0.2 gives q>0.5013.
    if b < a / 16.0 || b + 0.5 * a >= 0.2 {
        return f64::NAN;
    }
    let sc = sinhc_half(a);
    let d = rawmul(Pair::from(a), sc);
    let q = rawmul(
        Pair::new(
            f64::from_bits(0x40040d931ff62706),
            -f64::from_bits(0x3caa6a0d6f814637),
        ),
        rawsum(Pair::from(b), Pair::new(0.5 * d.hi, 0.5 * d.lo)),
    );
    if !(q.hi >= 1e-12 && q.hi <= 0.5) {
        return f64::NAN;
    }
    let r = rawdiv(d, q);
    let z = rawmul(r, r);
    if !(z.hi <= 0.5) {
        return f64::NAN;
    }
    let u = rawdiv(rawmul(q, normal_deferred(z)), sc);
    let w = rawmul(u, u);
    let a2 = rawmul(Pair::from(a), Pair::from(a));
    let [p0, p1] = conversion_pair(&conversion_0, &conversion_1, w.hi);
    let [p2, p3] = conversion_pair(&conversion_2, &conversion_3, w.hi);
    let [p4, p5] = conversion_pair(&conversion_4, &conversion_5, w.hi);
    let [p6, p7] = conversion_pair(&conversion_6, &conversion_7, w.hi);
    let [p8, p9] = conversion_pair(&conversion_8, &conversion_9, w.hi);
    let rows = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9];
    let mut y = rows[9];
    for i in (0..9).rev() {
        y = fma(y, a2.hi, rows[i]);
    }
    let c = rawsum(Pair::from(1.0), rawmul(w, Pair::from(y)));
    renorm(rawmul(u, c)).value()
}

#[inline]
fn atan_rel_squared(q: f64) -> f64 {
    // Rounded degree-five Remez coefficients. The exact rational approximation
    // and source-order rounding bounds are checked by verify_as1_atan.py.
    let mut p = f64::from_bits(0xbfb6c9831b2f66d6);
    for c in [
        f64::from_bits(0x3fbc709e301f27e7),
        f64::from_bits(0xbfc24923ee7916eb),
        f64::from_bits(0x3fc999999947d463),
        f64::from_bits(0xbfd5555555554dde),
        1.0,
    ] {
        p = fma(p, q, c);
    }
    p
}
#[inline]
fn slope(h: f64, t: f64) -> f64 {
    let x = h - t;
    let rule: &[[[f64; 2]; 2]] = if x >= 20.0 {
        &RULE2
    } else if x >= 12.0 {
        &RULE1
    } else {
        &RULE3
    };
    let tt = t * t;
    let kk = fma(h, h, -tt);
    let T = 4.0 * tt;
    let mut terms = [0.0; 4];
    for (i, atom) in rule.iter().enumerate() {
        let z = atom[0][0];
        let w = atom[1][0];
        let A = kk + z;
        let inv = 1.0 / A;
        let q = (T * z) * (inv * inv);
        terms[i % 4] = fma(w * inv, atan_rel_squared(q), terms[i % 4]);
    }
    (terms[0] + terms[1]) + (terms[2] + terms[3])
}
fn frexp(x: f64) -> (f64, i32) {
    let e = ilogb(x) + 1;
    (scalbn(x, -e), e)
}
#[inline]
fn logratio(a: f64, b: f64) -> Pair {
    let (ma, ea) = frexp(a);
    let (mb, eb) = frexp(b);
    let r = ma / mb;
    let e = fma(-r, mb, ma) / (r * mb);
    let k = (ea - eb) as f64;
    let c = f64::from_bits(0x3fe62e42fefa39ef);
    let cl = f64::from_bits(0x3c7abc9e3b39803f);
    let p = k * c;
    rawsum(Pair::new(p, fma(k, c, -p) + k * cl), Pair::new(r.ln(), e))
}
pub(super) fn balanced(a: f64, b: f64) -> f64 {
    let lr = rawsum(
        logratio(a, b),
        Pair::new(
            -f64::from_bits(0x3fed67f1c864beb5),
            f64::from_bits(0x3c865b5a1b7ff5df),
        ),
    );
    let L = 2.0 * lr.hi;
    if !(L > 70.0) {
        return f64::NAN;
    }
    let l = L.ln();
    let r = 1.0 / L;
    let A = a * a / 4.0;
    let p1 = fma(9.0, l, -6.0 - A);
    let y0 = fma(p1, r, fma(-3.0, l, L));
    let c2 = fma(fma(13.5, l, -fma(3.0, A, 45.0)), l, fma(5.0, A, 39.0));
    let y = fma(c2, r * r, y0);
    let h = y.sqrt();
    let t = a / (2.0 * h);
    if !(h - t >= 7.0 && a >= 1e-100 && h <= 40.0 && t <= 0.7) {
        return f64::NAN;
    }
    let D = slope(h, t);
    let Al = fma(a / 2.0, a / 2.0, -A);
    let t2 = A / y;
    let t2l = (fma(-t2, y, A) + Al) / y;
    let g = rawsum(
        rawsum(rawsum(lr, Pair::from((D / h).ln())), Pair::from(-0.5 * y)),
        Pair::new(-0.5 * t2, -0.5 * t2l),
    );
    let gv = g.value();
    let p = -0.5 / (y * D);
    let l = -0.5 + (t2 - 3.0) / (2.0 * y);
    let lp = (1.5 - t2) / (y * y);
    let H = l - p;
    let H3 = fma(H, l - 2.0 * p, lp);
    let n = -gv / p;
    let dy = n * (1.0 + 0.5 * H * n) / (1.0 + n * (H + H3 * n / 6.0));
    let Y = rawsum(Pair::from(y), Pair::from(dy));
    let hh = Y.hi.sqrt();
    let hl = (fma(-hh, hh, Y.hi) + Y.lo) / (2.0 * hh);
    let s = a / hh;
    let sl = (fma(-s, hh, a) - s * hl) / hh;
    s + sl
}

#[inline]
fn atm_rational(p: &[Pair], d: &[Pair], x: Pair, last: usize) -> Pair {
    // Evaluate padded numerator and denominator independently in two lanes.
    let n = p.len().max(d.len());
    let get = |i: usize| {
        let p = p.get(i).copied().unwrap_or(Pair::from(0.0));
        let d = d.get(i).copied().unwrap_or(Pair::from(0.0));
        (F64x2::new(p.hi, d.hi), F64x2::new(p.lo, d.lo))
    };
    let h = F64x2::splat(x.hi);
    let l = F64x2::splat(x.lo);
    let mut v = get(n - 1).0;
    for i in (last..n - 1).rev() {
        v = v.mul_add(h, get(i).0);
    }
    let mut vl = F64x2::splat(0.0);
    for i in (0..last).rev() {
        let (c, cl) = get(i);
        let prod = v * h;
        let ep = v.mul_add(h, -prod);
        let s = prod + c;
        let z = s - prod;
        let es = (prod - (s - z)) + (c - z);
        vl = vl.mul_add(h, v.mul_add(l, (ep + es) + cl));
        v = s;
    }
    let h = v.to_array();
    let l = vl.to_array();
    let inv = 1.0 / h[1];
    let q = h[0] * inv;
    let er = fma(-q, h[1], h[0]) + l[0];
    Pair::new(q, fma(-q, l[1], er) * inv)
}
pub(super) fn atm(b: f64) -> f64 {
    if !(b > 0.0 && b < 1.0) {
        return f64::NAN;
    }
    if b < 2.0 * 0.3413447460685429 {
        let q = 0.5 * b;
        let sq = q * q;
        let h = f64::from_bits(0x3fbdd4020da64ff1) - sq;
        let low = ((f64::from_bits(0x3fbdd4020da64ff1) - h) - sq) + fma(-q, q, sq);
        let r = atm_rational(&shift_central_p, &shift_central_d, Pair::new(h, low), 2);
        return fma(b, r.hi, b * r.lo);
    }
    let p = fma(-0.5, b, 0.5);
    let q = -p.ln();
    let rh = q.sqrt();
    let rl = fma(-rh, rh, q) / (2.0 * rh);
    let v = if rh < 2.05 {
        atm_rational(&shift_tail1_p, &shift_tail1_d, Pair::new(rh - 1.75, rl), 4)
    } else if rh < 3.41 {
        atm_rational(&shift_tail2_p, &shift_tail2_d, Pair::new(rh - 2.75, rl), 2)
    } else {
        atm_rational(&shift_tail3_p, &shift_tail3_d, Pair::new(rh - 4.75, rl), 2)
    };
    -2.0 * v.value()
}
