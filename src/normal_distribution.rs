use crate::erf_cody::erfc_cody;

const NORM_CDF_ASYMPTOTIC_EXPANSION_FIRST_THRESHOLD: f64 = -10.0;
const NORM_CDF_ASYMPTOTIC_EXPANSION_SECOND_THRESHOLD: f64 = -67108864.0;
// 1.0 / f64::sqrt(f64::EPSILON);
const FRAC_SQRT_2_PI: f64 = 0.3989422804014327;

pub(crate) fn norm_pdf(x: f64) -> f64 {
    FRAC_SQRT_2_PI * (-0.5 * x * x).exp()
}

pub(crate) fn norm_cdf(z: f64) -> f64 {
    if z <= NORM_CDF_ASYMPTOTIC_EXPANSION_FIRST_THRESHOLD {
        let mut sum = 1.0;
        if z >= NORM_CDF_ASYMPTOTIC_EXPANSION_SECOND_THRESHOLD {
            let zsqr = z * z;
            let mut i = 1.0;
            let mut g = 1.0;
            let mut x;
            let mut y;
            let mut a = f64::MAX;
            let mut lasta;
            loop {
                lasta = a;
                x = (4.0 * i - 3.0) / zsqr;
                y = x * ((4.0 * i - 1.0) / zsqr);
                a = g * (x - y);
                sum -= a;
                g *= y;
                i += 1.0;
                a = a.abs();
                if !(lasta > a && a >= (sum * f64::EPSILON).abs()) {
                    break;
                }
            }
        }
        return -norm_pdf(z) * sum / z;
    }
    0.5 * erfc_cody(-z * (1.0 / std::f64::consts::SQRT_2))
}


pub(crate) fn inverse_norm_cdf(u: f64) -> f64 {
    //
    // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
    //
    // Produces the normal deviate Z corresponding to a given lower
    // tail area of u; Z is accurate to about 1 part in 10**16.
    // see http://lib.stat.cmu.edu/apstat/241
    //
    let split1 = 0.425;
    let split2 = 5.0;
    let const1 = 0.180625;
    let const2 = 1.6;

    // Coefficients for P close to 0.5
    let a = [3.387_132_872_796_366_5, 1.331_416_678_917_843_8E2, 1.971_590_950_306_551_3E3,
        1.373_169_376_550_946E4, 4.592_195_393_154_987E4, 6.726_577_092_700_87E4,
        3.343_057_558_358_813E4, 2.509_080_928_730_122_7E3];
    let b = [4.231_333_070_160_091E1, 6.871_870_074_920_579E2,
        5.394_196_021_424_751E3, 2.121_379_430_158_659_7E4,
        3.930_789_580_009_271E4, 2.872_908_573_572_194_3E4,
        5.226_495_278_852_854E3];
    // Coefficients for P not close to 0, 0.5 or 1.
    let c = [1.423_437_110_749_683_5, 4.630_337_846_156_546,
        5.769_497_221_460_691, 3.647_848_324_763_204_5,
        1.270_458_252_452_368_4, 2.417_807_251_774_506E-1,
        2.272_384_498_926_918_4E-2, 7.745_450_142_783_414E-4];
    let d = [2.053_191_626_637_759, 1.676_384_830_183_803_8,
        6.897_673_349_851E-1, 1.481_039_764_274_800_8E-1,
        1.519_866_656_361_645_7E-2, 5.475_938_084_995_345E-4,
        1.050_750_071_644_416_9E-9];
    // Coefficients for P very close to 0 or 1
    let e = [6.657_904_643_501_103, 5.463_784_911_164_114,
        1.784_826_539_917_291_3, 2.965_605_718_285_048_7E-1,
        2.653_218_952_657_612_4E-2, 1.242_660_947_388_078_4E-3,
        2.711_555_568_743_487_6E-5, 2.010_334_399_292_288_1E-7];
    let f = [5.998_322_065_558_88E-1, 1.369_298_809_227_358E-1,
        1.487_536_129_085_061_5E-2, 7.868_691_311_456_133E-4,
        1.846_318_317_510_054_8E-5, 1.421_511_758_316_446E-7,
        2.044_263_103_389_939_7E-15];

    if u <= 0.0 {
        return u.ln();
    } else if u >= 1.0 {
        return (1.0 - u).ln();
    }

    let q = u - 0.5;
    if q.abs() <= split1 {
        let r = const1 - q * q;
        let mut numerator = a[0];
        for a in a.iter().skip(1) {
            numerator = numerator * r + a;
        }
        let mut denominator = b[0];
        for b in b.iter().skip(1) {
            denominator = denominator * r + b;
        }
        q * numerator / (denominator * r + 1.0)
    } else {
        let mut r = if q < 0.0 { u } else { 1.0 - u };
        r = r.ln().abs().sqrt();
        let ret;
        if r <= split2 {
            r -= const2;
            let mut numerator = c[0];
            for c in c.iter().skip(1) {
                numerator = numerator * r + c;
            }
            let mut denominator = d[0];
            for d in d.iter().skip(1) {
                denominator = denominator * r + d;
            }
            ret = numerator / (denominator * r + 1.0);
        } else {
            r -= split2;
            let mut numerator = e[0];
            for e in e.iter().skip(1) {
                numerator = numerator * r + e;
            }
            let mut denominator = f[0];
            for f in f.iter().skip(1) {
                denominator = denominator * r + f;
            }
            ret = numerator / (denominator * r + 1.0);
        }
        if q < 0.0 { -ret } else { ret }
    }
}

