//! experimental large-moneyness route (10 < a < 1417), with two independent SIMD lanes.
//! Compensated rounding points are preserved; inverse-normal seeds use fused Horner.
use crate::experimental::lanes::F64x2;
use crate::experimental::large_constants::{ASYM, CAP_TABLE, EXP_C, LN2, LOGK, M_COEFF};
use crate::experimental::math::inverse_log_probability;

#[derive(Clone, Copy, Default)]
pub(super) struct DD {
    pub(super) hi: f64,
    pub(super) lo: f64,
}

#[inline(always)]
fn dd(hi: f64) -> DD {
    DD { hi, lo: 0.0 }
}
#[inline(always)]
fn two(a: f64, b: f64) -> DD {
    let s = a + b;
    let v = s - a;
    DD {
        hi: s,
        lo: (a - (s - v)) + (b - v),
    }
}
#[inline(always)]
fn neg(a: DD) -> DD {
    DD {
        hi: -a.hi,
        lo: -a.lo,
    }
}
#[inline(always)]
fn add(a: DD, b: DD) -> DD {
    let t = two(a.hi, b.hi);
    let e = (a.lo + b.lo) + t.lo;
    two(t.hi, e)
}
#[inline(always)]
fn sub(a: DD, b: DD) -> DD {
    add(a, neg(b))
}
#[inline(always)]
fn mul(a: DD, b: DD) -> DD {
    let h = a.hi * b.hi;
    let mut l = a.hi.mul_add(b.hi, -h);
    l = a.hi.mul_add(b.lo, l);
    l = a.lo.mul_add(b.hi, l);
    l = a.lo.mul_add(b.lo, l);
    two(h, l)
}
#[inline(always)]
fn scale(a: DD, b: f64) -> DD {
    mul(a, dd(b))
}
#[inline(always)]
fn divide(a: DD, b: DD) -> DD {
    let q = a.hi / b.hi;
    let r = sub(a, scale(b, q));
    two(q, (r.hi + r.lo) / b.hi)
}
#[inline(always)]
fn sqrtdd(a: DD) -> DD {
    let q = a.hi.sqrt();
    let r = sub(a, mul(dd(q), dd(q)));
    two(q, (r.hi + r.lo) / (2.0 * q))
}
#[inline(always)]
fn val(a: DD) -> f64 {
    a.hi + a.lo
}

// Positive normal or subnormal input, the same exact mantissa as frexp.
#[inline]
fn frexp(x: f64) -> (f64, i32) {
    if x == 0.0 || !x.is_finite() {
        return (x, 0);
    }
    let mut v = x;
    let mut correction = 0;
    if v.abs() < f64::MIN_POSITIVE {
        v *= 18014398509481984.0;
        correction = -54;
    }
    let bits = v.to_bits();
    let e = ((bits >> 52) & 2047) as i32 - 1022 + correction;
    (
        f64::from_bits((bits & 0x800fffffffffffff) | (1022u64 << 52)),
        e,
    )
}

#[inline]
fn scalbn(x: f64, n: i32) -> f64 {
    // Calls below only scale a normal price upward by at most 1021, or a
    // small nonnegative integer downward by fewer than 249 binary powers.
    debug_assert!((-1022..=1023).contains(&n));
    x * f64::from_bits(((n + 1023) as u64) << 52)
}

#[inline]
fn logdd(a: DD) -> DD {
    let (mut m, mut e) = frexp(a.hi);
    if m < f64::from_bits(0x3fe6a09e667f3bcd) {
        m *= 2.0;
        e -= 1;
    }
    let t = two(m.ln(), a.lo / a.hi);
    add(scale(LN2, e as f64), t)
}

// The dispatcher has already checked that the price is positive and normal.
// Its DD low part is zero, so logdd(dd(b)) uses this exact frexp mantissa and
// exponent with a positive zero correction. Keep its compensated sum unchanged.
#[inline]
fn log_positive_normal_price(b: f64) -> DD {
    debug_assert!(b.is_normal() && b > 0.0);
    let bits = b.to_bits();
    let mut e = ((bits >> 52) & 2047) as i32 - 1022;
    let mut m = f64::from_bits((bits & 0x800fffffffffffff) | (1022u64 << 52));
    if m < f64::from_bits(0x3fe6a09e667f3bcd) {
        m *= 2.0;
        e -= 1;
    }
    let t = two(m.ln(), 0.0);
    add(scale(LN2, e as f64), t)
}

#[inline(always)]
fn raw_stage(p: DD, z: DD, c: DD) -> DD {
    let ph = p.hi * z.hi;
    let pl = p.hi.mul_add(z.hi, -ph);
    let s = two(ph, c.hi);
    let mut l = (pl + s.lo) + c.lo;
    l = p.lo.mul_add(z.hi, l);
    l = p.hi.mul_add(z.lo, l);
    DD { hi: s.hi, lo: l }
}

#[inline]
fn reduced_cap(r: DD) -> DD {
    let j = 32.0_f64.mul_add(r.hi, 0.5) as usize;
    let c = j as f64 * 0.03125;
    let z = two(r.hi - c, r.lo);
    let mut h = EXP_C[12].hi;
    for k in (6..12).rev() {
        h = h.mul_add(z.hi, EXP_C[k].hi);
    }
    let mut p = dd(h);
    for k in (0..6).rev() {
        p = raw_stage(p, z, EXP_C[k]);
    }
    p = two(p.hi, p.lo);
    mul(p, CAP_TABLE[j - 16])
}

#[inline]
fn mills_fast(z: f64) -> f64 {
    if z < 16.0 {
        let k = (2.0 * z).floor() as i32;
        if k < -12 {
            return f64::NAN;
        }
        let c = &M_COEFF[(k + 12) as usize];
        let x = z - (2 * k + 1) as f64 * 0.25;
        let mut p = c[24].hi;
        for i in (0..24).rev() {
            p = p.mul_add(x, c[i].hi);
        }
        p
    } else {
        let q = 1.0 / (z * z);
        let mut p = ASYM[13];
        for i in (0..13).rev() {
            p = p.mul_add(q, ASYM[i]);
        }
        p / z
    }
}

#[inline(always)]
fn fast_pair_impl<const A: bool, const B: bool>(a: f64, b: f64) -> [f64; 2] {
    let ka = if A { (2.0 * a).floor() as i32 } else { 0 };
    let kb = if B { (2.0 * b).floor() as i32 } else { 0 };
    let ca = &M_COEFF[(ka + 12) as usize];
    let cb = &M_COEFF[(kb + 12) as usize];
    let x = F64x2::new(
        if A {
            a - (2 * ka + 1) as f64 * 0.25
        } else {
            1.0 / (a * a)
        },
        if B {
            b - (2 * kb + 1) as f64 * 0.25
        } else {
            1.0 / (b * b)
        },
    );
    let coef = |i: usize| {
        F64x2::new(
            if A {
                ca[i].hi
            } else if i <= 13 {
                ASYM[i]
            } else {
                0.0
            },
            if B {
                cb[i].hi
            } else if i <= 13 {
                ASYM[i]
            } else {
                0.0
            },
        )
    };
    let n = if A || B { 24 } else { 13 };
    let mut p = coef(n);
    for i in (0..n).rev() {
        p = p.mul_add(x, coef(i));
    }
    if !(A && B) {
        p = p / F64x2::new(if A { 1.0 } else { a }, if B { 1.0 } else { b });
    }
    p.to_array()
}

#[inline]
fn mills_fast_pair(a: f64, b: f64) -> [f64; 2] {
    if a < -6.0 || b < -6.0 {
        return [f64::NAN; 2];
    }
    match (a < 16.0, b < 16.0) {
        (true, true) => fast_pair_impl::<true, true>(a, b),
        (true, false) => fast_pair_impl::<true, false>(a, b),
        (false, true) => fast_pair_impl::<false, true>(a, b),
        (false, false) => fast_pair_impl::<false, false>(a, b),
    }
}

#[inline(always)]
fn precise_pair_impl<const A: bool, const B: bool>(a: DD, b: DD) -> [DD; 2] {
    // Pack independent lanes only; retain every rounding point of the source.
    let prepare = |mut z: DD, tab: bool| {
        let k = if tab { (2.0 * z.hi).floor() as i32 } else { 0 };
        let mut iv = dd(1.0);
        let mut q = dd(0.0);
        if !(A && B) {
            let zh = if tab { 1.0 } else { z.hi };
            let ih = 1.0 / zh;
            let il = (-zh).mul_add(ih, 1.0) * ih;
            iv = DD { hi: ih, lo: il };
            let qh = ih * ih;
            q = DD {
                hi: qh,
                lo: ih.mul_add(ih, -qh) + (ih + ih) * il,
            };
        }
        if tab {
            let local = two(z.hi, -(2 * k + 1) as f64 * 0.25);
            z.lo += local.lo;
            q = dd(local.hi);
        }
        (z, q, iv, &M_COEFF[(k + 12) as usize])
    };
    let (za, qa, ia, ca) = prepare(a, A);
    let (zb, qb, ib, cb) = prepare(b, B);
    let qh = F64x2::new(qa.hi, qb.hi);
    let ql = F64x2::new(qa.lo, qb.lo);
    let coef_hi = |i: usize| {
        F64x2::new(
            if A {
                ca[i].hi
            } else if i <= 18 {
                ASYM[i]
            } else {
                0.0
            },
            if B {
                cb[i].hi
            } else if i <= 18 {
                ASYM[i]
            } else {
                0.0
            },
        )
    };
    let coef_lo = |i: usize| {
        F64x2::new(
            if A { ca[i].lo } else { 0.0 },
            if B { cb[i].lo } else { 0.0 },
        )
    };
    let n = if A || B { 24 } else { 18 };
    let mut h = coef_hi(n);
    for i in (3..n).rev() {
        h = h.mul_add(qh, coef_hi(i));
    }
    let mut l = F64x2::splat(0.0);
    for i in (0..3).rev() {
        let ph = h * qh;
        let pl = h.mul_add(qh, -ph);
        let c = coef_hi(i);
        let sh = ph + c;
        let v = sh - ph;
        let sl = (ph - (sh - v)) + (c - v);
        let mut next = (pl + sl) + coef_lo(i);
        next = l.mul_add(qh, next);
        next = h.mul_add(ql, next);
        h = sh;
        l = next;
    }
    if !(A && B) {
        let ih = F64x2::new(ia.hi, ib.hi);
        let il = F64x2::new(ia.lo, ib.lo);
        let ph = h * ih;
        let mut pl = h.mul_add(ih, -ph);
        pl = h.mul_add(il, pl);
        pl = l.mul_add(ih, pl);
        h = ph;
        l = pl;
    }
    let der = F64x2::new(za.hi, zb.hi).mul_add(h, F64x2::splat(-1.0));
    let corr = der * F64x2::new(za.lo, zb.lo);
    let tail = l + corr;
    let sh = h + tail;
    let v = sh - h;
    let sl = (h - (sh - v)) + (tail - v);
    let [ha, hb] = sh.to_array();
    let [la, lb] = sl.to_array();
    [DD { hi: ha, lo: la }, DD { hi: hb, lo: lb }]
}

#[inline]
fn mills_precise_pair(a: DD, b: DD) -> [DD; 2] {
    if a.hi < -1.0 || b.hi < -1.0 {
        return [DD {
            hi: f64::NAN,
            lo: f64::NAN,
        }; 2];
    }
    match (a.hi < 16.0, b.hi < 16.0) {
        (true, true) => precise_pair_impl::<true, true>(a, b),
        (true, false) => precise_pair_impl::<true, false>(a, b),
        (false, true) => precise_pair_impl::<false, true>(a, b),
        (false, false) => precise_pair_impl::<false, false>(a, b),
    }
}

#[derive(Clone, Copy, Default)]
struct Precise {
    residual: DD,
    y: DD,
    s: DD,
    derivative: f64,
}

#[inline]
fn geometry(a: f64, x: f64) -> (DD, DD, f64, f64) {
    let xh = x * x;
    let xl = x.mul_add(x, -xh);
    let w = two(xh, 2.0 * a);
    let wl = w.lo + xl;
    let q = w.hi.sqrt();
    let ql = ((-q).mul_add(q, w.hi) + wl) / (2.0 * q);
    let y = two(q, ql);
    let s = if x >= 0.0 {
        let t = two(y.hi, x);
        two(t.hi, t.lo + y.lo)
    } else {
        let den = two(y.hi, -x);
        let dl = den.lo + y.lo;
        let inv = 1.0 / den.hi;
        let sh = (2.0 * a) * inv;
        let sl = (-sh).mul_add(dl, (-sh).mul_add(den.hi, 2.0 * a)) * inv;
        two(sh, sl)
    };
    (y, s, xh, xl)
}

#[inline]
fn residual(xh: f64, xl: f64, d: DD, target: DD) -> DD {
    let bits = d.hi.to_bits();
    let mut e = ((bits >> 52) & 2047) as i32 - 1022;
    let mut m = f64::from_bits((bits & 0x000fffffffffffff) | (1022u64 << 52));
    if m < f64::from_bits(0x3fe6a09e667f3bcd) {
        m *= 2.0;
        e -= 1;
    }
    let lm = m.ln();
    let eh = e as f64 * LN2.hi;
    let el = (e as f64).mul_add(LN2.hi, -eh) + e as f64 * LN2.lo;
    let mut t = two(-0.5 * xh, eh);
    let mut low = (-0.5 * xl + el) + t.lo;
    t = two(t.hi, lm);
    low += t.lo;
    t = two(t.hi, -LOGK.hi);
    low += t.lo - LOGK.lo;
    t = two(t.hi, -target.hi);
    low += t.lo - target.lo;
    low += d.lo / d.hi;
    two(t.hi, low)
}

#[inline]
fn precise(a: f64, x: f64, target: DD, upper: bool) -> Precise {
    let (y, s, xh, xl) = geometry(a, x);
    let mm = mills_precise_pair(dd(if upper { x } else { -x }), y);
    let d = if upper {
        add(mm[0], mm[1])
    } else {
        sub(mm[0], mm[1])
    };
    let r = residual(xh, xl, d, target);
    let derivative = (if upper { -1.0 } else { 1.0 }) * val(s) / val(y) / val(d);
    Precise {
        residual: r,
        y,
        s,
        derivative,
    }
}

struct Eval {
    f: f64,
    df: f64,
    h: f64,
    y: f64,
}

#[inline]
fn eval_fast(a: f64, x: f64, target: f64, upper: bool) -> Eval {
    let y = x.mul_add(x, 2.0 * a).sqrt();
    let s = if x < 0.0 { 2.0 * a / (y - x) } else { x + y };
    let mm = mills_fast_pair(if upper { x } else { -x }, y);
    let d = if upper { mm[0] + mm[1] } else { mm[0] - mm[1] };
    let f = (-0.5 * x * x - LOGK.hi + d.ln()) - target;
    let df = (if upper { -1.0 } else { 1.0 }) * (s / y) / d;
    Eval {
        f,
        df,
        h: -x + 1.0 / y - x / (y * y),
        y,
    }
}

/// experimental defaults: early precise handoff, log seed, scalar finish and transport.
pub(crate) fn solve(a: f64, b: f64) -> f64 {
    if !(a > 10.0 && a < 1417.0 && b.is_normal() && b > 0.0 && a.is_finite()) {
        return f64::NAN;
    }
    let ell = add(log_positive_normal_price(b), dd(a / 2.0));
    let upper = ell.hi > -0.6;
    let mut target = ell;
    if upper {
        if ell.hi <= -f64::from_bits(0x3eb0000000000000) {
            let e = ell.hi.exp_m1();
            target = logdd(two(-e, -e.mul_add(ell.lo, ell.lo)));
        } else {
            let k = (a / (2.0 * LN2.hi)) as i32 - 1;
            let bs = scalbn(b, k);
            if !(bs > 0.0 && bs < 1.0 && bs.is_finite()) {
                return f64::NAN;
            }
            let r = sub(dd(a / 2.0), scale(LN2, k as f64));
            let mut cap = reduced_cap(r);
            let mut gap = sub(cap, dd(bs));
            if gap.hi.abs() <= f64::from_bits(0x3af0000000000000) {
                match fixed_scaled_cap(a, bs, k) {
                    Some((g, c)) => {
                        gap = g;
                        cap = c;
                    }
                    None => return f64::NAN,
                }
            }
            if !(gap.hi > 0.0) {
                return f64::NAN;
            }
            target = logdd(divide(gap, cap));
        }
        if !(target.hi < -0.7 && target.hi > -100.0) {
            return f64::NAN;
        }
    }
    let q = inverse_log_probability(val(target));
    let mut x = if upper { -q } else { q };
    let mut left = if upper { 0.0 } else { -40.0 };
    let mut right = if upper { 16.0 } else { 1.0 };
    if !(x > left && x < right) {
        x = (left + right) / 2.0;
    }
    let mut steps = 0;
    if upper {
        let yy = x.mul_add(x, 2.0 * a).sqrt();
        let ss = x + yy;
        let next = x + physical_step7(x, yy, ss);
        if next.is_finite() && next > left && next < right {
            x = next;
        }
        steps = 1;
    }
    let mut pe = Precise::default();
    let mut converged = false;
    while steps < 208 {
        if steps >= 1 {
            pe = precise(a, x, target, upper);
            if pe.residual.hi.is_finite()
                && pe.residual.hi.abs() <= f64::from_bits(0x3e90000000000000)
            {
                converged = true;
                break;
            }
        }
        let e = eval_fast(a, x, val(target), upper);
        if e.f.is_finite() && e.f.abs() <= f64::from_bits(0x3dd0000000000000) {
            pe = precise(a, x, target, upper);
            converged = true;
            break;
        }
        if !e.f.is_finite() {
            return f64::NAN;
        }
        if (e.f < 0.0) != upper {
            left = x;
        } else {
            right = x;
        }
        let mut nx = x + log_step5(a, x, &e);
        if steps % 4 == 3 || !(nx > left && nx < right) || !nx.is_finite() {
            nx = (left + right) / 2.0;
        }
        x = nx;
        steps += 1;
    }
    if !converged {
        return f64::NAN;
    }
    let yy = val(pe.y);
    let d = pe.derivative;
    let iy = 1.0 / yy;
    let iysq = iy * iy;
    let h = -x + iy - x * iysq;
    let hp = -1.0 - x * iy * iysq + (x * x - 2.0 * a) * (iysq * iysq);
    let aa = h - d;
    let bb = hp + aa * (h - 2.0 * d);
    let n = -val(pe.residual) / d;
    let numerator = (0.5 * aa).mul_add(n, 1.0);
    let denominator = n.mul_add((bb / 6.0).mul_add(n, aa), 1.0);
    let final_step = dd(n * (numerator / denominator));
    let corrected = add(dd(x), final_step);
    let dx = val(final_step);
    if dx.is_finite() && dx.abs() <= f64::from_bits(0x3ed0000000000000) {
        let invy = 1.0 / val(pe.y);
        let z = dx * invy;
        let k = x * invy;
        let c = 0.5 * (1.0 - k);
        let inc = z.mul_add((c * z).mul_add((-k).mul_add(z, 1.0), 1.0), 0.0);
        let low = pe.s.hi.mul_add(inc, pe.s.lo);
        return pe.s.hi + low;
    }
    let yf = sqrtdd(add(mul(corrected, corrected), dd(2.0 * a)));
    let sf = if corrected.hi < 0.0 {
        divide(dd(2.0 * a), sub(yf, corrected))
    } else {
        add(yf, corrected)
    };
    val(sf)
}

// Integer interval refinement on the original 2^-248 fixed grid.
#[derive(Clone, Copy, Default, PartialEq, Eq)]
struct U256([u64; 4]);

impl U256 {
    fn one() -> Self {
        Self([0, 0, 0, 1u64 << 56])
    }
    fn small(x: u64) -> Self {
        Self([x, 0, 0, 0])
    }
    fn cmp(self, b: Self) -> std::cmp::Ordering {
        for i in (0..4).rev() {
            if self.0[i] != b.0[i] {
                return self.0[i].cmp(&b.0[i]);
            }
        }
        std::cmp::Ordering::Equal
    }
    fn add(self, b: Self) -> Self {
        let mut r = Self::default();
        let mut carry = 0u128;
        for i in 0..4 {
            carry = self.0[i] as u128 + b.0[i] as u128 + (carry >> 64);
            r.0[i] = carry as u64;
        }
        debug_assert_eq!(carry >> 64, 0);
        r
    }
    fn sub(self, b: Self) -> Self {
        debug_assert!(self.cmp(b) != std::cmp::Ordering::Less);
        let mut r = Self::default();
        let mut borrow = 0u128;
        for i in 0..4 {
            let q = (1u128 << 64) + self.0[i] as u128 - b.0[i] as u128 - borrow;
            r.0[i] = q as u64;
            borrow = 1 - (q >> 64);
        }
        debug_assert_eq!(borrow, 0);
        r
    }
    fn shr(self, n: u32) -> Self {
        let mut r = Self::default();
        if n >= 256 {
            return r;
        }
        let k = (n / 64) as usize;
        let s = n % 64;
        for i in 0..4 - k {
            r.0[i] = self.0[i + k] >> s;
            if s != 0 && i + k + 1 < 4 {
                r.0[i] |= self.0[i + k + 1] << (64 - s);
            }
        }
        r
    }
    fn low_nonzero(self, n: u32) -> bool {
        if n >= 256 {
            return self != Self::default();
        }
        for i in 0..(n / 64) as usize {
            if self.0[i] != 0 {
                return true;
            }
        }
        let s = n % 64;
        s != 0 && (self.0[(n / 64) as usize] & ((1u64 << s) - 1)) != 0
    }
    fn ceil_shr(self, n: u32) -> Self {
        let r = self.shr(n);
        if self.low_nonzero(n) {
            r.add(Self::small(1))
        } else {
            r
        }
    }
    fn bits(self) -> u32 {
        for i in (0..4).rev() {
            if self.0[i] != 0 {
                return 64 * i as u32 + 64 - self.0[i].leading_zeros();
            }
        }
        0
    }
    fn mul_grid(self, b: Self) -> Self {
        let mut p = [0u64; 8];
        for i in 0..4 {
            let mut carry = 0u128;
            for j in 0..4 {
                let z = self.0[i] as u128 * b.0[j] as u128 + p[i + j] as u128 + carry;
                p[i + j] = z as u64;
                carry = z >> 64;
            }
            p[i + 4] = carry as u64;
        }
        debug_assert_eq!(p[7] >> 56, 0);
        let mut r = Self::default();
        for i in 0..4 {
            r.0[i] = (p[i + 3] >> 56) | (p[i + 4] << 8);
        }
        r
    }
    fn div_small(self, d: u32) -> Self {
        let mut r = Self::default();
        let mut rem = 0u128;
        for i in (0..4).rev() {
            let z = (rem << 64) | self.0[i] as u128;
            r.0[i] = (z / d as u128) as u64;
            rem = z % d as u128;
        }
        r
    }
    fn exact_double(x: f64) -> Self {
        debug_assert!(x > 0.0 && x.is_normal());
        let bits = x.to_bits();
        let e = ((bits >> 52) & 2047) as i32 - 1023;
        let mant = (bits & ((1u64 << 52) - 1)) | (1u64 << 52);
        let shift = e - 52 + 248;
        debug_assert!(shift >= 0 && shift + 53 <= 256);
        let mut r = Self::default();
        let k = (shift / 64) as usize;
        let s = shift % 64;
        r.0[k] = mant << s;
        if s != 0 && k + 1 < 4 {
            r.0[k + 1] = mant >> (64 - s);
        }
        r
    }
    fn to_double(self) -> f64 {
        let n = self.bits();
        if n == 0 {
            return 0.0;
        }
        if n <= 53 {
            return scalbn(self.0[0] as f64, -248);
        }
        let shift = n - 53;
        let top = self.shr(shift);
        let mut mant = top.0[0];
        let half = ((self.0[((shift - 1) / 64) as usize] >> ((shift - 1) % 64)) & 1) != 0;
        if half && (self.low_nonzero(shift - 1) || (mant & 1) != 0) {
            mant += 1;
        }
        scalbn(mant as f64, shift as i32 - 248)
    }
}

fn reduced_grid(a: f64, k: u32) -> U256 {
    const LN2_Q: [u64; 4] = [
        0x2d8a0d175b8baafa,
        0xaf40f343267298b6,
        0xabc9e3b39803f2f6,
        0x00b17217f7d1cf79,
    ];
    let mut ah = [0u64; 5];
    let mut lh = [0u64; 5];
    let mut rr = [0u64; 5];
    let bits = a.to_bits();
    let e = ((bits >> 52) & 2047) as i32 - 1023;
    let mant = (bits & ((1u64 << 52) - 1)) | (1u64 << 52);
    let shift = (e - 52 + 248 - 1) as u32;
    let w = (shift / 64) as usize;
    let t = shift % 64;
    debug_assert!(shift + 53 <= 320);
    ah[w] = mant << t;
    if t != 0 {
        ah[w + 1] = mant >> (64 - t);
    }
    let mut carry = 0u128;
    for i in 0..4 {
        let v = LN2_Q[i] as u128 * k as u128 + carry;
        lh[i] = v as u64;
        carry = v >> 64;
    }
    lh[4] = carry as u64;
    let mut borrow = 0u128;
    for i in 0..5 {
        let v = (1u128 << 64) + ah[i] as u128 - lh[i] as u128 - borrow;
        rr[i] = v as u64;
        borrow = 1 - (v >> 64);
    }
    debug_assert!(borrow == 0 && rr[4] == 0);
    U256([rr[0], rr[1], rr[2], rr[3]])
}

fn grid_dd(x: U256) -> DD {
    let h = x.to_double();
    if h == 0.0 {
        return DD::default();
    }
    let u = U256::exact_double(h);
    let l = if x.cmp(u) != std::cmp::Ordering::Less {
        x.sub(u).to_double()
    } else {
        -u.sub(x).to_double()
    };
    DD { hi: h, lo: l }
}

fn fixed_scaled_cap(a: f64, bs: f64, k: i32) -> Option<(DD, DD)> {
    let y = reduced_grid(a, k as u32).shr(5);
    let mut t = U256::one();
    let mut c = U256::one();
    for j in 1..=48 {
        t = t.mul_grid(y).div_small(j);
        // Zero is absorbing in the integer recurrence; all later terms are zero.
        if t == U256::default() {
            break;
        }
        c = if (j & 1) != 0 { c.sub(t) } else { c.add(t) };
    }
    for _ in 0..5 {
        c = c.mul_grid(c);
    }
    let radius = U256::small(8192);
    let cl = c.sub(radius);
    let ch = c.add(radius);
    let b = U256::exact_double(bs);
    if b.cmp(ch) != std::cmp::Ordering::Less || b.cmp(cl) != std::cmp::Ordering::Less {
        return None;
    }
    Some((grid_dd(c.sub(b)), grid_dd(c)))
}

/// The shared v40 fixed-cap route. Status values preserve the source enum:
/// positive=0, negative=1, unresolved=2, outside_window=3, bad_input=4.
pub(crate) fn resolve_cap(a: f64, b: f64) -> (u32, f64) {
    use std::cmp::Ordering::{Greater, Less};
    if !(a >= f64::from_bits(0x3d70000000000000)
        && a <= 10.0
        && b >= 0.00390625
        && b <= 1.0
        && b.is_finite())
    {
        return (4, f64::NAN);
    }
    let y = U256::exact_double(a * 0.03125);
    let mut t = U256::one();
    let mut c = U256::one();
    for k in 1..=48 {
        t = t.mul_grid(y).div_small(k);
        // Zero is absorbing in the integer recurrence; all later terms are zero.
        if t == U256::default() {
            break;
        }
        c = if (k & 1) != 0 { c.sub(t) } else { c.add(t) };
    }
    for _ in 0..4 {
        c = c.mul_grid(c);
    }
    let radius = U256::small(2048);
    let cl = c.sub(radius);
    let ch = c.add(radius);
    let p = U256::exact_double(b);
    if ch.cmp(p) == Less {
        return (1, f64::NAN);
    }
    if cl.cmp(p) != Greater {
        return (2, f64::NAN);
    }
    let g = c.sub(p);
    let gl = g.sub(radius);
    let gh = g.add(radius);
    if gl.cmp(ch.ceil_shr(201)) == Less {
        return (2, f64::NAN);
    }
    if gh.cmp(cl.shr(87)) == Greater {
        return (3, f64::NAN);
    }
    (0, g.to_double() / c.to_double())
}

#[inline(always)]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    a.mul_add(b, c)
}

#[inline]
#[allow(non_snake_case, unused_variables)]
fn physical_step7(x: f64, y: f64, s: f64) -> f64 {
    let n = mills_fast(y) * (y / s);
    let iy = f64::from_bits(0x3ff0000000000000) / y;
    let k = x * iy;
    let v = n * iy;
    let nn = n * n;
    let xn = -x * n;
    let g0 = f64::from_bits(0x3ff0000000000000);
    let g1 = xn;
    let w0 = f64::from_bits(0x3ff0000000000000);
    let g2 = fma(xn, g1, -nn * g0) / f64::from_bits(0x4000000000000000);
    let g3 = fma(xn, g2, -nn * g1) / f64::from_bits(0x4008000000000000);
    let g4 = fma(xn, g3, -nn * g2) / f64::from_bits(0x4010000000000000);
    let g5 = fma(xn, g4, -nn * g3) / f64::from_bits(0x4014000000000000);
    let g6 = fma(xn, g5, -nn * g4) / f64::from_bits(0x4018000000000000);
    let v1 = v;
    let w1 = (f64::from_bits(0x3ff0000000000000) - k) * (f64::from_bits(0x3ff0000000000000)) * v1;
    let v2 = v1 * v / f64::from_bits(0x4000000000000000);
    let w2 = (f64::from_bits(0x3ff0000000000000) - k)
        * (fma(
            -f64::from_bits(0x4008000000000000),
            k,
            f64::from_bits(0x0000000000000000),
        ))
        * v2;
    let v3 = v2 * v / f64::from_bits(0x4008000000000000);
    let w3 = (f64::from_bits(0x3ff0000000000000) - k)
        * (fma(
            fma(
                f64::from_bits(0x402e000000000000),
                k,
                f64::from_bits(0x0000000000000000),
            ),
            k,
            -f64::from_bits(0x4008000000000000),
        ))
        * v3;
    let v4 = v3 * v / f64::from_bits(0x4010000000000000);
    let w4 = (f64::from_bits(0x3ff0000000000000) - k)
        * (fma(
            fma(
                fma(
                    -f64::from_bits(0x405a400000000000),
                    k,
                    f64::from_bits(0x0000000000000000),
                ),
                k,
                f64::from_bits(0x4046800000000000),
            ),
            k,
            f64::from_bits(0x0000000000000000),
        ))
        * v4;
    let v5 = v4 * v / f64::from_bits(0x4014000000000000);
    let w5 = (f64::from_bits(0x3ff0000000000000) - k)
        * (fma(
            fma(
                fma(
                    fma(
                        f64::from_bits(0x408d880000000000),
                        k,
                        f64::from_bits(0x0000000000000000),
                    ),
                    k,
                    -f64::from_bits(0x4083b00000000000),
                ),
                k,
                f64::from_bits(0x0000000000000000),
            ),
            k,
            f64::from_bits(0x4046800000000000),
        ))
        * v5;
    let v6 = v5 * v / f64::from_bits(0x4018000000000000);
    let w6 = (f64::from_bits(0x3ff0000000000000) - k)
        * (fma(
            fma(
                fma(
                    fma(
                        fma(
                            -f64::from_bits(0x40c44d8000000000),
                            k,
                            f64::from_bits(0x0000000000000000),
                        ),
                        k,
                        f64::from_bits(0x40c2750000000000),
                    ),
                    k,
                    f64::from_bits(0x0000000000000000),
                ),
                k,
                -f64::from_bits(0x40989c0000000000),
            ),
            k,
            f64::from_bits(0x0000000000000000),
        ))
        * v6;
    let q2 = (fma(g0, w1, g1)) / f64::from_bits(0x4000000000000000);
    let q3 = (fma(g0, w2, fma(g1, w1, g2))) / f64::from_bits(0x4008000000000000);
    let q4 = (fma(g0, w3, fma(g1, w2, fma(g2, w1, g3)))) / f64::from_bits(0x4010000000000000);
    let q5 = (fma(g0, w4, fma(g1, w3, fma(g2, w2, fma(g3, w1, g4)))))
        / f64::from_bits(0x4014000000000000);
    let q6 = (fma(
        g0,
        w5,
        fma(g1, w4, fma(g2, w3, fma(g3, w2, fma(g4, w1, g5)))),
    )) / f64::from_bits(0x4018000000000000);
    let q7 = (fma(
        g0,
        w6,
        fma(
            g1,
            w5,
            fma(g2, w4, fma(g3, w3, fma(g4, w2, fma(g5, w1, g6)))),
        ),
    )) / f64::from_bits(0x401c000000000000);
    let c0 = f64::from_bits(0x3ff0000000000000);
    let c1 = f64::from_bits(0x3ff0000000000000);
    let c2 = fma(q2, c0, c1);
    let c3 = fma(q3, c0, fma(q2, c1, c2));
    let c4 = fma(q4, c0, fma(q3, c1, fma(q2, c2, c3)));
    let c5 = fma(q5, c0, fma(q4, c1, fma(q3, c2, fma(q2, c3, c4))));
    let c6 = fma(
        q6,
        c0,
        fma(q5, c1, fma(q4, c2, fma(q3, c3, fma(q2, c4, c5)))),
    );
    let c7 = fma(
        q7,
        c0,
        fma(
            q6,
            c1,
            fma(q5, c2, fma(q4, c3, fma(q3, c4, fma(q2, c5, c6)))),
        ),
    );
    n * c6 / c7
}

#[inline]
#[allow(non_snake_case, unused_variables)]
fn log_step5(a: f64, x: f64, e: &Eval) -> f64 {
    let n = -e.f / e.df;
    let d = e.df;
    let H = e.h;
    let iy = f64::from_bits(0x3ff0000000000000) / e.y;
    let iy2 = iy * iy;
    let xx = x * x * iy2;
    let H1 = -f64::from_bits(0x3ff0000000000000) - x * iy * iy2
        + (x * x - f64::from_bits(0x4000000000000000) * a) * iy2 * iy2;
    let H2 = (-f64::from_bits(0x3ff0000000000000) + f64::from_bits(0x4008000000000000) * xx)
        * iy
        * iy2
        + (f64::from_bits(0x4018000000000000) * x - f64::from_bits(0x4020000000000000) * x * xx)
            * iy2
            * iy2;
    let H3 = (f64::from_bits(0x4048000000000000) * xx * xx
        - f64::from_bits(0x4048000000000000) * xx
        + f64::from_bits(0x4018000000000000))
        * iy2
        * iy2
        + x * (f64::from_bits(0x4022000000000000) - f64::from_bits(0x402e000000000000) * xx)
            * iy
            * iy2
            * iy2;
    let A = H - d;
    let B = H1 + A * (H - f64::from_bits(0x4000000000000000) * d);
    let C = H2
        + (f64::from_bits(0x4008000000000000) * H - f64::from_bits(0x4010000000000000) * d) * H1
        + A * (H * H - f64::from_bits(0x4018000000000000) * H * d
            + f64::from_bits(0x4018000000000000) * d * d);
    let D = H3
        + (f64::from_bits(0x4010000000000000) * H - f64::from_bits(0x4014000000000000) * d) * H2
        + f64::from_bits(0x4008000000000000) * H1 * H1
        + H1 * (f64::from_bits(0x4034000000000000) * d * d
            - f64::from_bits(0x4039000000000000) * d * H
            + f64::from_bits(0x4018000000000000) * H * H)
        + (((H - f64::from_bits(0x402e000000000000) * d) * H
            + f64::from_bits(0x4049000000000000) * d * d)
            * H
            - f64::from_bits(0x404e000000000000) * d * d * d)
            * H
        + f64::from_bits(0x4038000000000000) * d * d * d * d;
    let num = fma(
        n,
        fma(
            n,
            fma(
                C / f64::from_bits(0x4038000000000000),
                n,
                A * A / f64::from_bits(0x4010000000000000) + B / f64::from_bits(0x4008000000000000),
            ),
            f64::from_bits(0x3ff8000000000000) * A,
        ),
        f64::from_bits(0x3ff0000000000000),
    );
    let den = fma(
        n,
        fma(
            n,
            fma(
                n,
                fma(
                    D / f64::from_bits(0x405e000000000000),
                    n,
                    A * B / f64::from_bits(0x4018000000000000)
                        + C / f64::from_bits(0x4028000000000000),
                ),
                f64::from_bits(0x4008000000000000) * A * A / f64::from_bits(0x4010000000000000)
                    + B / f64::from_bits(0x4000000000000000),
            ),
            f64::from_bits(0x4000000000000000) * A,
        ),
        f64::from_bits(0x3ff0000000000000),
    );
    n * num / den
}

#[cfg(test)]
mod tests {
    use super::{dd, log_positive_normal_price, logdd};

    #[test]
    fn normal_price_log_preserves_both_compensated_components() {
        // Sweep normal exponents and the mantissa-renormalization seam.
        // The general DD logarithm remains the reference for gap/ratio inputs.
        let cutoff = 0x0006_a09e_667f_3bcd_u64;
        let fractions = [
            0,
            1,
            cutoff - 1,
            cutoff,
            cutoff + 1,
            (1 << 52) - 2,
            (1 << 52) - 1,
        ];
        for exponent in 1_u64..2047 {
            for fraction in fractions {
                let b = f64::from_bits((exponent << 52) | fraction);
                let expected = logdd(dd(b));
                let actual = log_positive_normal_price(b);
                assert_eq!(actual.hi.to_bits(), expected.hi.to_bits(), "b={b:e}");
                assert_eq!(actual.lo.to_bits(), expected.lo.to_bits(), "b={b:e}");
            }
        }
    }
}
