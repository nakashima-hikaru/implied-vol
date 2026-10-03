//! Independent mathematical-root fixtures for the source-preserving experimental backend.
//! High-precision core roots: supplied archive `input/all_references.json`.
//! Tiny/ATM/seam roots: original repository `fresh_accuracy.json` and
//! `optimization/fresh_combined_fma.json` (direct Black, 140/440 digits).
//! Paper roots: archived `reference_{dataset}.bin` hi+lo values, not solver outputs.
//! These are finite regression checks of rho <= 1, not a global precision proof.
#![cfg(feature = "experimental")]

use implied_vol::solver::{BlackSolver, Experimental};
use implied_vol::{ImpliedBlackVolatility, ImpliedBlackVolatilityNormalised};

struct Reference {
    name: &'static str,
    a: f64,
    b: f64,
    root_hi: f64,
    root_lo: f64,
    eta: f64,
}

const REFERENCES: &[Reference] = &[
    // a=0x1.4000000000001p+3, b=0x1.0000000000000p-1022
    Reference {
        name: "just above a=10, minimum normal price",
        a: f64::from_bits(0x4024_0000_0000_0001),
        b: f64::from_bits(0x0010_0000_0000_0000),
        root_hi: f64::from_bits(0x3fd1_1e3c_aff7_356b),
        root_lo: f64::from_bits(0x3c60_0cbc_831d_5065),
        eta: f64::from_bits(0x3cb0_02ec_8fc3_a097),
    },
    // a=0x1.4000000000000p+5, b=0x1.0000000000000p-1022
    Reference {
        name: "a=40, minimum normal price",
        a: f64::from_bits(0x4044_0000_0000_0000),
        b: f64::from_bits(0x0010_0000_0000_0000),
        root_hi: f64::from_bits(0x3ff1_1a52_de0d_04c9),
        root_lo: f64::from_bits(0xbc8a_0fc1_39d5_6692),
        eta: f64::from_bits(0x3cb0_02eb_5eca_357a),
    },
    // a=0x1.2c00000000000p+8, b=0x1.0000000000000p-1022
    Reference {
        name: "a=300, minimum normal price",
        a: f64::from_bits(0x4072_c000_0000_0000),
        b: f64::from_bits(0x0010_0000_0000_0000),
        root_hi: f64::from_bits(0x4020_1a24_d92e_29ff),
        root_lo: f64::from_bits(0xbccf_6bbf_163c_1431),
        eta: f64::from_bits(0x3cb0_02fa_6e48_1422),
    },
    // a=0x1.6200000000000p+10, b=0x1.0000000000000p-1022
    Reference {
        name: "a=1416, minimum normal price",
        a: f64::from_bits(0x4096_2000_0000_0000),
        b: f64::from_bits(0x0010_0000_0000_0000),
        root_hi: f64::from_bits(0x404a_d7a5_5d16_706a),
        root_lo: f64::from_bits(0xbce3_7a13_bae8_7029),
        eta: f64::from_bits(0x3cb0_8f6c_9c34_5fa6),
    },
    // a=0x1.6232bdd7abcd2p+10, b=0x1.0000000000000p-1022
    Reference {
        name: "highest supported sampled moneyness",
        a: f64::from_bits(0x4096_232b_dd7a_bcd2),
        b: f64::from_bits(0x0010_0000_0000_0000),
        root_hi: f64::from_bits(0x404e_a651_1e0f_1f29),
        root_lo: f64::from_bits(0xbcec_5c4f_ee83_e192),
        eta: f64::from_bits(0x3ef4_9ae8_19df_55ce),
    },
    // a=0x1.e8bdbfcd9144ep+4, b=0x1.f3e558cf4de54p-23
    Reference {
        name: "fixed-grid cap refinement with extraordinarily small gap",
        a: f64::from_bits(0x403e_8bdb_fcd9_144e),
        b: f64::from_bits(0x3e8f_3e55_8cf4_de54),
        root_hi: f64::from_bits(0x403a_545b_24df_27a0),
        root_lo: f64::from_bits(0x3cdc_b1de_895b_f01f),
        eta: f64::from_bits(0x42f7_9a4a_3443_e12e),
    },
    // a=0x1.4000000000001p+3, b=0x1.b993fe00d536fp-8
    Reference {
        name: "cap near a=10",
        a: f64::from_bits(0x4024_0000_0000_0001),
        b: f64::from_bits(0x3f7b_993f_e00d_536f),
        root_hi: f64::from_bits(0x4032_0362_c944_7aad),
        root_lo: f64::from_bits(0xbc8b_10f5_e1a2_4f8f),
        eta: f64::from_bits(0x3fb9_9c35_1401_b182),
    },
    // a=0x1.2c00000000000p+8, b=0x1.82e16284f5ec4p-217
    Reference {
        name: "cap near a=300",
        a: f64::from_bits(0x4072_c000_0000_0000),
        b: f64::from_bits(0x3268_2e16_284f_5ec4),
        root_hi: f64::from_bits(0x4041_0f54_e4c2_edae),
        root_lo: f64::from_bits(0xbc9f_0ec9_d32d_5b9b),
        eta: f64::from_bits(0x3f87_34e8_f86c_f045),
    },
    // a=0x1.6200000000000p+10, b=0x1.7c8ab2288c9abp-1022
    Reference {
        name: "cap near a=1416",
        a: f64::from_bits(0x4096_2000_0000_0000),
        b: f64::from_bits(0x0017_c8ab_2288_c9ab),
        root_hi: f64::from_bits(0x404f_172c_00e7_960a),
        root_lo: f64::from_bits(0xbcd6_a7d4_e637_f2ae),
        eta: f64::from_bits(0x3f83_6e6a_0d45_5ccd),
    },
    // a=0x1.df3273cea4c32p+3, b=0x1.424d551f7a9a0p-12
    Reference {
        name: "a~15, log relative price~-0.599999",
        a: f64::from_bits(0x402d_f327_3cea_4c32),
        b: f64::from_bits(0x3f34_24d5_51f7_a9a0),
        root_hi: f64::from_bits(0x4017_1b39_2bab_f813),
        root_lo: f64::from_bits(0x3c9d_54a1_27ea_13b7),
        eta: f64::from_bits(0x3cb3_fb1b_f064_ab07),
    },
    // a=0x1.0000000000000p+4, b=0x1.7cd79b5647c95p-15
    Reference {
        name: "a~15, log relative price~-2.0",
        a: f64::from_bits(0x4030_0000_0000_0000),
        b: f64::from_bits(0x3f07_cd79_b564_7c95),
        root_hi: f64::from_bits(0x4013_4204_abd9_2660),
        root_lo: f64::from_bits(0xbca7_d696_0499_097a),
        eta: f64::from_bits(0x3cb1_b714_2731_18ac),
    },
    // a=0x1.246dbb6037facp+4, b=0x1.fc57d173fe998p-43
    Reference {
        name: "a~15, log relative price~-20.0",
        a: f64::from_bits(0x4032_46db_b603_7fac),
        b: f64::from_bits(0x3d4f_c57d_173f_e998),
        root_hi: f64::from_bits(0x4004_ef03_5c2c_4586),
        root_lo: f64::from_bits(0x3c90_cc54_6334_183b),
        eta: f64::from_bits(0x3cb0_51c7_f3f2_bf73),
    },
    // a=0x1.e027c8574aaefp+4, b=0x1.66c76614911eep-23
    Reference {
        name: "a~30, log relative price~-0.599999",
        a: f64::from_bits(0x403e_027c_8574_aaef),
        b: f64::from_bits(0x3e86_6c76_6149_11ee),
        root_hi: f64::from_bits(0x401f_ff28_e81e_ca0b),
        root_lo: f64::from_bits(0x3ca6_1fcd_3ef1_576d),
        eta: f64::from_bits(0x3cb2_d669_1f74_a590),
    },
    // a=0x1.0000000000000p+5, b=0x1.0576460718281p-26
    Reference {
        name: "a~30, log relative price~-2.0",
        a: f64::from_bits(0x4040_0000_0000_0000),
        b: f64::from_bits(0x3e50_5764_6071_8281),
        root_hi: f64::from_bits(0x401c_596d_c54e_28e1),
        root_lo: f64::from_bits(0xbbe2_37c0_f194_4686),
        eta: f64::from_bits(0x3cb1_3a0a_329a_3bc0),
    },
    // a=0x1.0a24f4acbd56ap+5, b=0x1.0b23d96c81bf0p-53
    Reference {
        name: "a~30, log relative price~-20.0",
        a: f64::from_bits(0x4040_a24f_4acb_d56a),
        b: f64::from_bits(0x3ca0_b23d_96c8_1bf0),
        root_hi: f64::from_bits(0x4010_ef1c_fc53_54ee),
        root_lo: f64::from_bits(0xbcb5_7046_15c1_c2c1),
        eta: f64::from_bits(0x3cb0_43c2_9af1_8e49),
    },
    // a=0x1.9000000000000p+6, b=0x1.ffde793b82cbbp-74
    Reference {
        name: "a~100, log relative price~-0.599999",
        a: f64::from_bits(0x4059_0000_0000_0000),
        b: f64::from_bits(0x3b5f_fde7_93b8_2cbb),
        root_hi: f64::from_bits(0x402c_ac0e_52c1_5b2c),
        root_lo: f64::from_bits(0x3cc0_16bb_4bac_2297),
        eta: f64::from_bits(0x3cb1_9067_835c_c8b3),
    },
    // a=0x1.9000000000000p+6, b=0x1.f8e6c24b5590fp-76
    Reference {
        name: "a~100, log relative price~-2.0",
        a: f64::from_bits(0x4059_0000_0000_0000),
        b: f64::from_bits(0x3b3f_8e6c_24b5_590f),
        root_hi: f64::from_bits(0x402a_4d40_37f5_ce6f),
        root_lo: f64::from_bits(0xbcca_c436_18d1_5d56),
        eta: f64::from_bits(0x3cb0_b354_fb17_ac1e),
    },
    // a=0x1.b50791d34ab69p+6, b=0x1.c2255942c1153p-108
    Reference {
        name: "a~100, log relative price~-20.0",
        a: f64::from_bits(0x405b_5079_1d34_ab69),
        b: f64::from_bits(0x393c_2255_942c_1153),
        root_hi: f64::from_bits(0x4024_3a4c_19cb_3375),
        root_lo: f64::from_bits(0xbccb_0e63_3b18_01bc),
        eta: f64::from_bits(0x3cb0_2b15_59e0_d565),
    },
    // a=0x1.2ce01aa4a5c1ep+9, b=0x1.0b69bd46ee148p-435
    Reference {
        name: "a~600, log relative price~-0.599999",
        a: f64::from_bits(0x4082_ce01_aa4a_5c1e),
        b: f64::from_bits(0x24c0_b69b_d46e_e148),
        root_hi: f64::from_bits(0x4041_6bef_bb13_f954),
        root_lo: f64::from_bits(0x3ce8_d10a_fb5a_292f),
        eta: f64::from_bits(0x3cb0_a394_a914_b9c5),
    },
    // a=0x1.429c4e590ecafp+9, b=0x1.a0e9b27c01190p-469
    Reference {
        name: "a~600, log relative price~-2.0",
        a: f64::from_bits(0x4084_29c4_e590_ecaf),
        b: f64::from_bits(0x22aa_0e9b_27c0_1190),
        root_hi: f64::from_bits(0x4041_6fba_39f8_0c91),
        root_lo: f64::from_bits(0x3cdd_1aaa_f1cd_e315),
        eta: f64::from_bits(0x3cb0_472b_2dbc_9d6a),
    },
    // a=0x1.e2e33a1aad91ap+8, b=0x1.df52134502760p-378
    Reference {
        name: "a~600, log relative price~-20.0",
        a: f64::from_bits(0x407e_2e33_a1aa_d91a),
        b: f64::from_bits(0x285d_f521_3450_2760),
        root_hi: f64::from_bits(0x4039_c923_17d3_2c92),
        root_lo: f64::from_bits(0xbcdc_5af5_5a18_454f),
        eta: f64::from_bits(0x3cb0_1578_7c4f_52e7),
    },
    // a=0x1.447bb0306cc85p+10, b=0x1.d525e83bda411p-938
    Reference {
        name: "a~1300, log relative price~-0.599999",
        a: f64::from_bits(0x4094_47bb_0306_cc85),
        b: f64::from_bits(0x055d_525e_83bd_a411),
        root_hi: f64::from_bits(0x4049_8bc9_7a82_3074),
        root_lo: f64::from_bits(0xbcc8_ceba_f0d3_9d5e),
        eta: f64::from_bits(0x3cb0_6f68_4637_e4c5),
    },
    // a=0x1.3e84bbfeee55ap+10, b=0x1.0c3322b990589p-922
    Reference {
        name: "a~1300, log relative price~-2.0",
        a: f64::from_bits(0x4093_e84b_bfee_e55a),
        b: f64::from_bits(0x0650_c332_2b99_0589),
        root_hi: f64::from_bits(0x4048_b487_71ba_fdc0),
        root_lo: f64::from_bits(0xbcb8_eb61_5613_eeb1),
        eta: f64::from_bits(0x3cb0_3280_0168_4c02),
    },
    // a=0x1.46b8106bebc10p+10, b=0x1.3e65ff95b8a1ep-972
    Reference {
        name: "a~1300, log relative price~-20.0",
        a: f64::from_bits(0x4094_6b81_06be_bc10),
        b: f64::from_bits(0x0333_e65f_f95b_8a1e),
        root_hi: f64::from_bits(0x4046_cb5b_0eb3_4cc0),
        root_lo: f64::from_bits(0xbcb7_f4d8_1fe6_9ebc),
        eta: f64::from_bits(0x3cb0_0d25_2382_95f3),
    },
    // a=0x1.5fdf5f2388d53p+4, b=0x1.d64c7d8fd35dfp-357
    Reference {
        name: "difficult near-unit volatility and deep lower tail",
        a: f64::from_bits(0x4035_fdf5_f238_8d53),
        b: f64::from_bits(0x29ad_64c7_d8fd_35df),
        root_hi: f64::from_bits(0x3ff0_12bd_d483_1963),
        root_lo: f64::from_bits(0xbc9f_ebab_3523_d01d),
        eta: f64::from_bits(0x3cb0_087f_8227_04f0),
    },
    Reference {
        name: "ATM generating total volatility 1.000000e-300",
        a: f64::from_bits(0x0000_0000_0000_0000),
        b: f64::from_bits(0x0191_194b_2f79_3d21),
        root_hi: f64::from_bits(0x01a5_6e1f_c2f8_f358),
        root_lo: f64::from_bits(0x0000_0000_00c2_27a9),
        eta: f64::from_bits(0x3cc0_0000_0000_0000),
    },
    Reference {
        name: "ATM generating total volatility 1.000000e-20",
        a: f64::from_bits(0x0000_0000_0000_0000),
        b: f64::from_bits(0x3bb2_d6ea_8e41_2acb),
        root_hi: f64::from_bits(0x3bc7_9ca1_0c92_4224),
        root_lo: f64::from_bits(0xb86a_23f4_ee8a_647e),
        eta: f64::from_bits(0x3cc0_0000_0000_0000),
    },
    Reference {
        name: "ATM generating total volatility 2.000000e-01",
        a: f64::from_bits(0x0000_0000_0000_0000),
        b: f64::from_bits(0x3fb4_6450_7526_870f),
        root_hi: f64::from_bits(0x3fc9_9999_9999_999a),
        root_lo: f64::from_bits(0xbc6f_648f_19f3_2068),
        eta: f64::from_bits(0x3cc0_06d7_207d_bfa7),
    },
    Reference {
        name: "ATM generating total volatility 8.000000e+00",
        a: f64::from_bits(0x0000_0000_0000_0000),
        b: f64::from_bits(0x3fef_ff7b_2943_55a8),
        root_hi: f64::from_bits(0x401f_ffff_ffff_ffff),
        root_lo: f64::from_bits(0x3c82_83f5_1f17_aef9),
        eta: f64::from_bits(0x3d4d_37ae_263e_2d83),
    },
    Reference {
        name: "microscopic a~1.000000e-300, s/a~1.0",
        a: f64::from_bits(0x01a5_6e1f_c2f8_f359),
        b: f64::from_bits(0x016c_9143_9e69_1be6),
        root_hi: f64::from_bits(0x01a5_6e1f_c2f8_f359),
        root_lo: f64::from_bits(0x8000_0000_0021_6e27),
        eta: f64::from_bits(0x3cb5_8256_2b0a_7af2),
    },
    Reference {
        name: "microscopic a~1.000000e-300, s/a~10.0",
        a: f64::from_bits(0x01a5_6e1f_c2f8_f359),
        b: f64::from_bits(0x01c2_cd2f_d9d7_fb0a),
        root_hi: f64::from_bits(0x01da_c9a7_b3b7_302f),
        root_lo: f64::from_bits(0x8000_0000_00ff_f60f),
        eta: f64::from_bits(0x3cbe_252a_86e9_d91b),
    },
    Reference {
        name: "microscopic a~1.000000e-110, s/a~1.0",
        a: f64::from_bits(0x2918_0c90_3f73_79f2),
        b: f64::from_bits(0x28e0_077e_f0b5_b2b3),
        root_hi: f64::from_bits(0x2918_0c90_3f73_79f2),
        root_lo: f64::from_bits(0x25a2_d50a_7c36_6c34),
        eta: f64::from_bits(0x3cb5_8256_2b0a_7af2),
    },
    Reference {
        name: "microscopic a~1.000000e-110, s/a~10.0",
        a: f64::from_bits(0x2918_0c90_3f73_79f2),
        b: f64::from_bits(0x2935_1963_9c05_f9f1),
        root_hi: f64::from_bits(0x294e_0fb4_4f50_586e),
        root_lo: f64::from_bits(0x25d7_a738_4893_ede4),
        eta: f64::from_bits(0x3cbe_252a_86e9_d91b),
    },
    Reference {
        name: "microscopic a~1.000000e-20, s/a~1.0",
        a: f64::from_bits(0x3bc7_9ca1_0c92_4223),
        b: f64::from_bits(0x3b8f_79c7_237e_6f8d),
        root_hi: f64::from_bits(0x3bc7_9ca1_0c92_4223),
        root_lo: f64::from_bits(0x3840_bad5_ad05_976b),
        eta: f64::from_bits(0x3cb5_8256_2b0a_7af2),
    },
    Reference {
        name: "microscopic a~1.000000e-20, s/a~10.0",
        a: f64::from_bits(0x3bc7_9ca1_0c92_4223),
        b: f64::from_bits(0x3be4_b72f_4e28_96fd),
        root_hi: f64::from_bits(0x3bfd_83c9_4fb6_d2ac),
        root_lo: f64::from_bits(0xb89d_4307_c6f5_5030),
        eta: f64::from_bits(0x3cbe_252a_86e9_d91b),
    },
    Reference {
        name: "a=0x1.a36e2eb1c432cp-14 inflection-side seam",
        a: f64::from_bits(0x3f1a_36e2_eb1c_432c),
        b: f64::from_bits(0x3fb4_7aec_6344_33d7),
        root_hi: f64::from_bits(0x3fc9_ba34_ade1_59e9),
        root_lo: f64::from_bits(0xbc6d_ab1b_a941_61e0),
        eta: f64::from_bits(0x3cc0_05a0_2ef4_54ac),
    },
    Reference {
        name: "a=0x1.a36e2eb1c432dp-14 inflection-side seam",
        a: f64::from_bits(0x3f1a_36e2_eb1c_432d),
        b: f64::from_bits(0x3fb4_7aec_6344_33d7),
        root_hi: f64::from_bits(0x3fc9_ba34_ade1_59e9),
        root_lo: f64::from_bits(0xbc6d_a108_e872_3625),
        eta: f64::from_bits(0x3cc0_05a0_2ef4_54ac),
    },
    Reference {
        name: "a=0x1.a36e2eb1c432ep-14 inflection-side seam",
        a: f64::from_bits(0x3f1a_36e2_eb1c_432e),
        b: f64::from_bits(0x3fb4_7aec_6344_33d7),
        root_hi: f64::from_bits(0x3fc9_ba34_ade1_59e9),
        root_lo: f64::from_bits(0xbc6d_96f6_27a3_0a69),
        eta: f64::from_bits(0x3cc0_05a0_2ef4_54ac),
    },
    Reference {
        name: "a=0x1.9999999999999p-4 inflection-side seam",
        a: f64::from_bits(0x3fb9_9999_9999_9999),
        b: f64::from_bits(0x3fc6_366f_0589_4c80),
        root_hi: f64::from_bits(0x3fe1_dd3e_fa77_4b9e),
        root_lo: f64::from_bits(0x3c69_46b6_4453_f5d0),
        eta: f64::from_bits(0x3cbd_2c07_d296_a972),
    },
    Reference {
        name: "a=0x1.999999999999ap-4 inflection-side seam",
        a: f64::from_bits(0x3fb9_9999_9999_999a),
        b: f64::from_bits(0x3fc6_366f_0589_4c80),
        root_hi: f64::from_bits(0x3fe1_dd3e_fa77_4b9e),
        root_lo: f64::from_bits(0x3c7e_b866_6e78_dda5),
        eta: f64::from_bits(0x3cbd_2c07_d296_a972),
    },
    Reference {
        name: "a=0x1.999999999999bp-4 inflection-side seam",
        a: f64::from_bits(0x3fb9_9999_9999_999b),
        b: f64::from_bits(0x3fc6_366f_0589_4c7f),
        root_hi: f64::from_bits(0x3fe1_dd3e_fa77_4b9e),
        root_lo: f64::from_bits(0xbc81_f8f0_2f6a_88a9),
        eta: f64::from_bits(0x3cbd_2c07_d296_a972),
    },
    Reference {
        name: "a=0x1.6666666666665p-2 inflection-side seam",
        a: f64::from_bits(0x3fd6_6666_6666_6665),
        b: f64::from_bits(0x3fcb_8d20_1564_42c5),
        root_hi: f64::from_bits(0x3fee_29e6_e289_91a2),
        root_lo: f64::from_bits(0xbc83_44ac_3752_e80a),
        eta: f64::from_bits(0x3cba_f6db_474c_c472),
    },
    Reference {
        name: "a=0x1.6666666666666p-2 inflection-side seam",
        a: f64::from_bits(0x3fd6_6666_6666_6666),
        b: f64::from_bits(0x3fcb_8d20_1564_42c5),
        root_hi: f64::from_bits(0x3fee_29e6_e289_91a2),
        root_lo: f64::from_bits(0x3c7b_d42f_27cf_2b1e),
        eta: f64::from_bits(0x3cba_f6db_474c_c472),
    },
    Reference {
        name: "a=0x1.6666666666667p-2 inflection-side seam",
        a: f64::from_bits(0x3fd6_6666_6666_6667),
        b: f64::from_bits(0x3fcb_8d20_1564_42c5),
        root_hi: f64::from_bits(0x3fee_29e6_e289_91a3),
        root_lo: f64::from_bits(0xbc80_e724_a0dd_ecd9),
        eta: f64::from_bits(0x3cba_f6db_474c_c472),
    },
    Reference {
        name: "a=0x1.fffffffffffffp-2 inflection-side seam",
        a: f64::from_bits(0x3fdf_ffff_ffff_ffff),
        b: f64::from_bits(0x3fcb_ef82_29d3_49c3),
        root_hi: f64::from_bits(0x3ff1_ae07_701c_218d),
        root_lo: f64::from_bits(0xbc86_3982_6c03_b1b4),
        eta: f64::from_bits(0x3cba_38e3_977b_8bef),
    },
    Reference {
        name: "a=0x1.0000000000000p-1 inflection-side seam",
        a: f64::from_bits(0x3fe0_0000_0000_0000),
        b: f64::from_bits(0x3fcb_ef82_29d3_49c3),
        root_hi: f64::from_bits(0x3ff1_ae07_701c_218d),
        root_lo: f64::from_bits(0x3c73_fed2_34c4_ab49),
        eta: f64::from_bits(0x3cba_38e3_977b_8bef),
    },
    Reference {
        name: "a=0x1.0000000000001p-1 inflection-side seam",
        a: f64::from_bits(0x3fe0_0000_0000_0001),
        b: f64::from_bits(0x3fcb_ef82_29d3_49c4),
        root_hi: f64::from_bits(0x3ff1_ae07_701c_218e),
        root_lo: f64::from_bits(0xbc4c_d631_ae38_d81f),
        eta: f64::from_bits(0x3cba_38e3_977b_8bef),
    },
    Reference {
        name: "a=0x1.fffffffffffffp-1 inflection-side seam",
        a: f64::from_bits(0x3fef_ffff_ffff_ffff),
        b: f64::from_bits(0x3fc9_6bd4_08e7_e758),
        root_hi: f64::from_bits(0x3ff8_48ae_a761_aadc),
        root_lo: f64::from_bits(0x3c3e_5498_7107_70de),
        eta: f64::from_bits(0x3cb8_b228_7f5e_6de9),
    },
    Reference {
        name: "a=0x1.0000000000000p+0 inflection-side seam",
        a: f64::from_bits(0x3ff0_0000_0000_0000),
        b: f64::from_bits(0x3fc9_6bd4_08e7_e757),
        root_hi: f64::from_bits(0x3ff8_48ae_a761_aadc),
        root_lo: f64::from_bits(0xbc63_7904_4049_1d25),
        eta: f64::from_bits(0x3cb8_b228_7f5e_6de8),
    },
    Reference {
        name: "a=0x1.0000000000001p+0 inflection-side seam",
        a: f64::from_bits(0x3ff0_0000_0000_0001),
        b: f64::from_bits(0x3fc9_6bd4_08e7_e757),
        root_hi: f64::from_bits(0x3ff8_48ae_a761_aadd),
        root_lo: f64::from_bits(0xbc77_2e07_169f_a321),
        eta: f64::from_bits(0x3cb8_b228_7f5e_6de8),
    },
    Reference {
        name: "a=0x1.7ffffffffffffp+1 inflection-side seam",
        a: f64::from_bits(0x4007_ffff_ffff_ffff),
        b: f64::from_bits(0x3fb6_ace7_45c4_58fe),
        root_hi: f64::from_bits(0x4004_6988_a190_f433),
        root_lo: f64::from_bits(0x3c8c_0ef3_8d94_b139),
        eta: f64::from_bits(0x3cb6_4560_e569_207d),
    },
    Reference {
        name: "a=0x1.8000000000000p+1 inflection-side seam",
        a: f64::from_bits(0x4008_0000_0000_0000),
        b: f64::from_bits(0x3fb6_ace7_45c4_58fb),
        root_hi: f64::from_bits(0x4004_6988_a190_f433),
        root_lo: f64::from_bits(0xbc86_c540_103c_9a08),
        eta: f64::from_bits(0x3cb6_4560_e569_207d),
    },
    Reference {
        name: "a=0x1.8000000000001p+1 inflection-side seam",
        a: f64::from_bits(0x4008_0000_0000_0001),
        b: f64::from_bits(0x3fb6_ace7_45c4_58fc),
        root_hi: f64::from_bits(0x4004_6988_a190_f434),
        root_lo: f64::from_bits(0x3c8f_b02f_4a43_dc45),
        eta: f64::from_bits(0x3cb6_4560_e569_207d),
    },
    Reference {
        name: "a=0x1.3ffffffffffffp+3 inflection-side seam",
        a: f64::from_bits(0x4023_ffff_ffff_ffff),
        b: f64::from_bits(0x3f69_1d26_1209_44cd),
        root_hi: f64::from_bits(0x4012_4b03_0e9b_aa7f),
        root_lo: f64::from_bits(0x3c73_0596_c6df_adfa),
        eta: f64::from_bits(0x3cb4_0294_01e7_8925),
    },
    Reference {
        name: "a=0x1.4000000000000p+3 inflection-side seam",
        a: f64::from_bits(0x4024_0000_0000_0000),
        b: f64::from_bits(0x3f69_1d26_1209_44ca),
        root_hi: f64::from_bits(0x4012_4b03_0e9b_aa80),
        root_lo: f64::from_bits(0x3c86_d17d_0336_7e3e),
        eta: f64::from_bits(0x3cb4_0294_01e7_8925),
    },
    Reference {
        name: "a=0x1.4000000000001p+3 inflection-side seam",
        a: f64::from_bits(0x4024_0000_0000_0001),
        b: f64::from_bits(0x3f69_1d26_1209_44c1),
        root_hi: f64::from_bits(0x4012_4b03_0e9b_aa80),
        root_lo: f64::from_bits(0xbc79_6485_8c48_bab2),
        eta: f64::from_bits(0x3cb4_0294_01e7_8925),
    },
    Reference {
        name: "microscopic smallest saved normal price",
        a: f64::from_bits(0x03d0_8004_0b48_08c4),
        b: f64::from_bits(0x0198_f99e_f8ce_b8a0),
        root_hi: f64::from_bits(0x03a5_e7cc_fe0b_253e),
        root_lo: f64::from_bits(0x0003_2b21_f3de_4378),
        eta: f64::from_bits(0x3cb0_6890_b2dd_1650),
    },
    Reference {
        name: "CLY-3D independent worst saved rho input #9769",
        a: f64::from_bits(0x3fe6_fed4_2880_3712),
        b: f64::from_bits(0x3fc0_a0e4_a08f_85d2),
        root_hi: f64::from_bits(0x3ff0_0e15_213e_b647),
        root_lo: f64::from_bits(0xbc71_a23b_6f77_4066),
        eta: f64::from_bits(0x3cb7_9c2e_cdea_48ca),
    },
    Reference {
        name: "CLY-20 independent worst saved rho input #330",
        a: f64::from_bits(0x3fc5_7df3_1504_0cf7),
        b: f64::from_bits(0x3f85_af0e_9df0_7fec),
        root_hi: f64::from_bits(0x3fc3_9dde_6e3c_7d92),
        root_lo: f64::from_bits(0xbc4b_8db4_58bb_f630),
        eta: f64::from_bits(0x3cb5_107a_0672_c1bc),
    },
    Reference {
        name: "CLY-80 independent worst saved rho input #477",
        a: f64::from_bits(0x3ff0_b81b_bfa8_6488),
        b: f64::from_bits(0x3fb7_ab82_92f4_7e76),
        root_hi: f64::from_bits(0x3ff1_a7cd_c69d_14d4),
        root_lo: f64::from_bits(0xbc85_749d_64a1_d127),
        eta: f64::from_bits(0x3cb6_2078_1fdc_917c),
    },
    Reference {
        name: "Jaeckel independent worst saved rho input #4836",
        a: f64::from_bits(0x4000_19ce_d3dc_5a67),
        b: f64::from_bits(0x3fd6_4238_e51b_207d),
        root_hi: f64::from_bits(0x4012_bee2_907f_eab4),
        root_lo: f64::from_bits(0xbcab_b794_174b_a213),
        eta: f64::from_bits(0x3cd0_b449_cf8a_ac13),
    },
    Reference {
        name: "Market independent worst saved rho input #2495",
        a: f64::from_bits(0x3fa0_59a8_14e8_11e4),
        b: f64::from_bits(0x3fc8_82f8_da9b_dc05),
        root_hi: f64::from_bits(0x3fe0_cccc_cccc_cccc),
        root_lo: f64::from_bits(0x3c74_f975_41e9_a35c),
        eta: f64::from_bits(0x3cbf_2b75_5e15_d53f),
    },
    Reference {
        name: "Stress independent worst saved rho input #782",
        a: f64::from_bits(0x3ffe_8d0e_c4ba_5a01),
        b: f64::from_bits(0x3fd1_fe33_8fcf_ccc7),
        root_hi: f64::from_bits(0x4009_0b94_c99f_da69),
        root_lo: f64::from_bits(0x3ca7_ffb4_af78_3b2a),
        eta: f64::from_bits(0x3cbe_c4c6_e58d_1675),
    },
    Reference {
        name: "HighVol independent worst saved rho input #25",
        a: f64::from_bits(0x4012_6bb1_bbb5_5515),
        b: f64::from_bits(0x3fa9_dbd8_ba46_618d),
        root_hi: f64::from_bits(0x400a_d533_6963_eefc),
        root_lo: f64::from_bits(0x3c91_15c9_0cdc_9663),
        eta: f64::from_bits(0x3cb6_5326_0791_06e8),
    },
];

fn inverse(x: f64, b: f64) -> Option<f64> {
    ImpliedBlackVolatilityNormalised::builder()
        .log_moneyness(x)
        .normalised_price(b)
        .build()
        .unwrap()
        .calculate_with::<Experimental>()
}

fn check_root(r: &Reference, actual: f64) {
    assert!(
        actual.is_finite() && actual > 0.0,
        "{}: invalid output {actual}",
        r.name
    );
    // Near the root, Sterbenz makes the leading subtraction exact. Keeping
    // the independent low part avoids replacing the mathematical root with
    // an already-rounded f64 target.
    let rho = (((actual - r.root_hi) - r.root_lo).abs() / r.root_hi) / r.eta;
    assert!(rho <= 1.0, "{}: attainable-accuracy ratio {rho}", r.name);
}

#[test]
fn independent_roots_span_seams_tiny_scales_and_cap_rounding() {
    for r in REFERENCES {
        let actual = inverse(-r.a, r.b).unwrap();
        check_root(r, actual);
        assert_eq!(
            actual.to_bits(),
            inverse(r.a, r.b).unwrap().to_bits(),
            "{}: reciprocal-strike symmetry",
            r.name
        );
    }
}

#[test]
fn experimental_owns_normalized_boundaries_and_atm_math() {
    for r in REFERENCES {
        let expected = Experimental::implied_total_volatility(-r.a, r.b).unwrap();
        let actual = inverse(-r.a, r.b).unwrap();
        assert_eq!(actual.to_bits(), expected.to_bits());
    }
}

#[test]
fn rounded_cap_is_classified_by_the_exact_experimental_boundary() {
    // exp(-1416/2) rounds downward: its binary64 value remains strictly
    // below the mathematical cap, and therefore has a finite implied root.
    let a = 1416.0_f64;
    let b = f64::from_bits(0x0017_c8ab_2288_c9ab);
    assert_eq!(b, (-0.5 * a).exp());
    let reference = REFERENCES.iter().find(|r| r.a == a && r.b == b).unwrap();
    check_root(reference, inverse(-a, b).unwrap());
    // exp(-1/2) rounds upward: this represented price is unattainable.
    let cap_above_exact = f64::from_bits(0x3fe3_68b2_fc6f_960a);
    assert_eq!(cap_above_exact, (-0.5_f64).exp());
    assert!(inverse(-1.0, cap_above_exact).is_none());
    assert!(inverse(-1.0, cap_above_exact.next_up()).is_none());
    assert_eq!(inverse(0.0, 1.0), Some(f64::INFINITY));
    assert_eq!(inverse(-1.0, 0.0), Some(0.0));
}

#[test]
fn full_api_preserves_atm_roots_and_call_put_time_value_normalization() {
    for r in REFERENCES.iter().filter(|r| r.a == 0.0) {
        for is_call in [false, true] {
            let actual = ImpliedBlackVolatility::builder()
                .forward(1.0)
                .strike(1.0)
                .expiry(1.0)
                .is_call(is_call)
                .option_price(r.b)
                .build()
                .unwrap()
                .calculate_with::<Experimental>()
                .unwrap();
            check_root(r, actual);
        }
    }
    // Dyadic prices retain the time value exactly when intrinsic value is
    // added. The same OTM leg must survive call/put parity and F/K exchange.
    for time_value in [0.0625_f64, 0.125, 0.5] {
        for expiry in [0.25_f64, 1.0, 4.0] {
            let mut results = Vec::new();
            for (forward, strike) in [(1.0_f64, 2.0_f64), (2.0, 1.0)] {
                for is_call in [false, true] {
                    let intrinsic = (if is_call {
                        forward - strike
                    } else {
                        strike - forward
                    })
                    .max(0.0);
                    let actual = ImpliedBlackVolatility::builder()
                        .forward(forward)
                        .strike(strike)
                        .expiry(expiry)
                        .is_call(is_call)
                        .option_price(intrinsic + time_value)
                        .build()
                        .unwrap()
                        .calculate_with::<Experimental>()
                        .unwrap();
                    assert!(actual.is_finite() && actual > 0.0);
                    results.push(actual.to_bits());
                }
            }
            assert!(
                results.iter().all(|r| *r == results[0]),
                "call/put parity or reciprocal F/K changed the result"
            );
        }
    }
}
