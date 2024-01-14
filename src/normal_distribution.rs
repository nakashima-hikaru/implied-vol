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
    const SPLIT1: f64 = 0.425;
    const SPLIT2: f64 = 5.0;
    const CONST1: f64 = 0.180625;
    const CONST2: f64 = 1.6;

    // Coefficients for P close to 0.5
    const A0: f64 = 3.387_132_872_796_366_5;
    const A1: f64 = 1.331_416_678_917_843_8E2;
    const A2: f64 = 1.971_590_950_306_551_3E3;
    const A3: f64 = 1.373_169_376_550_946E4;
    const A4: f64 = 4.592_195_393_154_987E4;
    const A5: f64 = 6.726_577_092_700_87E4;
    const A6: f64 = 3.343_057_558_358_813E4;
    const A7: f64 = 2.509_080_928_730_122_7E3;
    const B1: f64 = 4.231_333_070_160_091E1;
    const B2: f64 = 6.871_870_074_920_579E2;
    const B3: f64 = 5.394_196_021_424_751E3;
    const B4: f64 = 2.121_379_430_158_659_7E4;
    const B5: f64 = 3.930_789_580_009_271E4;
    const B6: f64 = 2.872_908_573_572_194_3E4;
    const B7: f64 = 5.226_495_278_852_854E3;
    // Coefficients for P not close to 0, 0.5 or 1.
    const C0: f64 = 1.423_437_110_749_683_5;
    const C1: f64 = 4.630_337_846_156_546;
    const C2: f64 = 5.769_497_221_460_691;
    const C3: f64 = 3.647_848_324_763_204_5;
    const C4: f64 = 1.270_458_252_452_368_4;
    const C5: f64 = 2.417_807_251_774_506E-1;
    const C6: f64 = 2.272_384_498_926_918_4E-2;
    const C7: f64 = 7.745_450_142_783_414E-4;
    const D1: f64 = 2.053_191_626_637_759;
    const D2: f64 = 1.676_384_830_183_803_8;
    const D3: f64 = 6.897_673_349_851E-1;
    const D4: f64 = 1.481_039_764_274_800_8E-1;
    const D5: f64 = 1.519_866_656_361_645_7E-2;
    const D6: f64 = 5.475_938_084_995_345E-4;
    const D7: f64 = 1.050_750_071_644_416_9E-9;
    // Coefficients for P very close to 0 or 1
    const E0: f64 = 6.657_904_643_501_103;
    const E1: f64 = 5.463_784_911_164_114;
    const E2: f64 = 1.784_826_539_917_291_3;
    const E3: f64 = 2.965_605_718_285_048_7E-1;
    const E4: f64 = 2.653_218_952_657_612_4E-2;
    const E5: f64 = 1.242_660_947_388_078_4E-3;
    const E6: f64 = 2.711_555_568_743_487_6E-5;
    const E7: f64 = 2.010_334_399_292_288_1E-7;
    const F1: f64 = 5.998_322_065_558_88E-1;
    const F2: f64 = 1.369_298_809_227_358E-1;
    const F3: f64 = 1.487_536_129_085_061_5E-2;
    const F4: f64 = 7.868_691_311_456_133E-4;
    const F5: f64 = 1.846_318_317_510_054_8E-5;
    const F6: f64 = 1.421_511_758_316_446E-7;
    const F7: f64 = 2.044_263_103_389_939_7E-15;

    if u <= 0.0 {
        return u.ln();
    } else if u >= 1.0 {
        return (1.0 - u).ln();
    }

    let q = u - 0.5;
    if q.abs() <= SPLIT1 {
        let r = CONST1 - q * q;
        q * (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0) /
            (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + 1.0)
    } else {
        let mut r = if q < 0.0 { u } else { 1.0 - u };
        r = (-r.ln()).sqrt();
        let ret =
        if r < SPLIT2 {
            r -= CONST2;
            (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0) /
                (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + 1.0)
        } else {
            r -= SPLIT2;
            (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0) /
                (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + 1.0)
        };
        if q < 0.0 { -ret } else { ret }
    }
}