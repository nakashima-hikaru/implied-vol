//! Independent exact-input roots for the reduced-degree Experimental polynomials.
//! AS1 references cover all three quadrature rules and tiny scaling. The Mills
//! divided-difference cases cover each reduced degree and adjacent prices.
//! Direct high-precision Black inversion at 120/200 or 440/700 decimal digits
//! was repeated for each exact binary64 input; the saved low root part retains
//! its displacement from the nearest binary64 value.
#![cfg(feature = "experimental")]

use implied_vol::solver::{BlackSolver, Experimental};

struct Reference {
    name: &'static str,
    a: f64,
    b: f64,
    root_hi: f64,
    root_lo: f64,
    eta: f64,
}

const REFERENCES: &[Reference] = &[
    Reference {
        name: "as1_rule_3_largest_geometry_ratio",
        a: f64::from_bits(0x4024_0000_0000_0000),
        b: f64::from_bits(0x3cd7_b816_42c1_163b),
        root_hi: f64::from_bits(0x3ff4_e7ee_6c77_9501),
        root_lo: f64::from_bits(0x3c81_3c1a_4211_fe06),
        eta: f64::from_bits(0x3cb0_430f_6873_6715),
    },
    Reference {
        name: "as1_rule_3_scaled",
        a: f64::from_bits(0x2918_0c90_3f73_79f2),
        b: f64::from_bits(0x21d8_5111_9c77_2e03),
        root_hi: f64::from_bits(0x28e0_0860_2a4c_fbf7),
        root_lo: f64::from_bits(0xa585_c107_3ea6_0714),
        eta: f64::from_bits(0x3cb0_1bdf_19b2_edb9),
    },
    Reference {
        name: "as1_rule_1_largest_geometry_ratio",
        a: f64::from_bits(0x4024_0000_0000_0000),
        b: f64::from_bits(0x386f_a99e_44db_7478),
        root_hi: f64::from_bits(0x3fe9_ccc9_9a00_cd86),
        root_lo: f64::from_bits(0xbc78_505a_2a88_bfbf),
        eta: f64::from_bits(0x3cb0_1a26_1bdc_ff8d),
    },
    Reference {
        name: "as1_rule_1_scaled",
        a: f64::from_bits(0x2918_0c90_3f73_79f2),
        b: f64::from_bits(0x162a_d89b_ac51_8cf6),
        root_hi: f64::from_bits(0x28d3_3d40_32c2_c7f5),
        root_lo: f64::from_bits(0xa569_eded_2564_4ace),
        eta: f64::from_bits(0x3cb0_0a2a_0550_181c),
    },
    Reference {
        name: "as1_rule_2_largest_geometry_ratio",
        a: f64::from_bits(0x4024_0000_0000_0000),
        b: f64::from_bits(0x2cc2_90ba_8978_8d3e),
        root_hi: f64::from_bits(0x3fdf_9c15_2bc3_0be5),
        root_lo: f64::from_bits(0xbc7b_80b2_c4a6_735d),
        eta: f64::from_bits(0x3cb0_09eb_c28a_df47),
    },
    Reference {
        name: "dpoly_degree_13_below",
        a: f64::from_bits(0x3fea_505b_162a_32c5),
        b: f64::from_bits(0x3d23_7120_a09c_155c),
        root_hi: f64::from_bits(0x3fbe_6025_51a8_c960),
        root_lo: f64::from_bits(0xbc18_0c6b_c02a_a4f7),
        eta: f64::from_bits(0x3cb0_5072_55d5_acef),
    },
    Reference {
        name: "dpoly_degree_13_exact",
        a: f64::from_bits(0x3fea_505b_162a_32c5),
        b: f64::from_bits(0x3d23_7120_a09c_155d),
        root_hi: f64::from_bits(0x3fbe_6025_51a8_c960),
        root_lo: f64::from_bits(0x3bfd_7e6a_be6d_0e1d),
        eta: f64::from_bits(0x3cb0_5072_55d5_acef),
    },
    Reference {
        name: "dpoly_degree_13_above",
        a: f64::from_bits(0x3fea_505b_162a_32c5),
        b: f64::from_bits(0x3d23_7120_a09c_155e),
        root_hi: f64::from_bits(0x3fbe_6025_51a8_c960),
        root_lo: f64::from_bits(0x3c23_65d0_8fb0_9602),
        eta: f64::from_bits(0x3cb0_5072_55d5_acef),
    },
    Reference {
        name: "dpoly_degree_14_below",
        a: f64::from_bits(0x3fdd_3ce9_f0a5_3937),
        b: f64::from_bits(0x3e22_fed5_44a5_e183),
        root_hi: f64::from_bits(0x3fb6_c005_8eea_c3a7),
        root_lo: f64::from_bits(0xbc2c_f098_5d37_1634),
        eta: f64::from_bits(0x3cb0_8c11_a9c1_c98e),
    },
    Reference {
        name: "dpoly_degree_14_exact",
        a: f64::from_bits(0x3fdd_3ce9_f0a5_3937),
        b: f64::from_bits(0x3e22_fed5_44a5_e184),
        root_hi: f64::from_bits(0x3fb6_c005_8eea_c3a7),
        root_lo: f64::from_bits(0xbc0f_e1bb_805d_6c50),
        eta: f64::from_bits(0x3cb0_8c11_a9c1_c98e),
    },
    Reference {
        name: "dpoly_degree_14_above",
        a: f64::from_bits(0x3fdd_3ce9_f0a5_3937),
        b: f64::from_bits(0x3e22_fed5_44a5_e185),
        root_hi: f64::from_bits(0x3fb6_c005_8eea_c3a7),
        root_lo: f64::from_bits(0x3c19_ff75_3a10_c016),
        eta: f64::from_bits(0x3cb0_8c11_a9c1_c98e),
    },
    Reference {
        name: "dpoly_degree_15_below",
        a: f64::from_bits(0x3fcf_5dcb_ff43_b156),
        b: f64::from_bits(0x3e78_d675_3e1c_afc7),
        root_hi: f64::from_bits(0x3fad_0715_7b2f_f649),
        root_lo: f64::from_bits(0xbbfc_73e6_2c3c_9b96),
        eta: f64::from_bits(0x3cb0_bf08_d90b_924f),
    },
    Reference {
        name: "dpoly_degree_15_exact",
        a: f64::from_bits(0x3fcf_5dcb_ff43_b156),
        b: f64::from_bits(0x3e78_d675_3e1c_afc8),
        root_hi: f64::from_bits(0x3fad_0715_7b2f_f649),
        root_lo: f64::from_bits(0x3c14_cb5e_36bd_b826),
        eta: f64::from_bits(0x3cb0_bf08_d90b_924f),
    },
    Reference {
        name: "dpoly_degree_15_above",
        a: f64::from_bits(0x3fcf_5dcb_ff43_b156),
        b: f64::from_bits(0x3e78_d675_3e1c_afc9),
        root_hi: f64::from_bits(0x3fad_0715_7b2f_f649),
        root_lo: f64::from_bits(0x3c28_59da_fc45_4b98),
        eta: f64::from_bits(0x3cb0_bf08_d90b_924f),
    },
    Reference {
        name: "dpoly_degree_16_below",
        a: f64::from_bits(0x3fc2_8084_6271_aa64),
        b: f64::from_bits(0x3f2e_b369_0113_5f19),
        root_hi: f64::from_bits(0x3fb0_2061_446f_fa9a),
        root_lo: f64::from_bits(0x3be6_2dae_8da6_8ddf),
        eta: f64::from_bits(0x3cb2_134c_d213_4a07),
    },
    Reference {
        name: "dpoly_degree_16_exact",
        a: f64::from_bits(0x3fc2_8084_6271_aa64),
        b: f64::from_bits(0x3f2e_b369_0113_5f1a),
        root_hi: f64::from_bits(0x3fb0_2061_446f_fa9a),
        root_lo: f64::from_bits(0x3c32_22b7_6ebe_f92a),
        eta: f64::from_bits(0x3cb2_134c_d213_4a07),
    },
    Reference {
        name: "dpoly_degree_16_above",
        a: f64::from_bits(0x3fc2_8084_6271_aa64),
        b: f64::from_bits(0x3f2e_b369_0113_5f1b),
        root_hi: f64::from_bits(0x3fb0_2061_446f_fa9a),
        root_lo: f64::from_bits(0x3c41_ca00_b488_5ef3),
        eta: f64::from_bits(0x3cb2_134c_d213_4a07),
    },
    Reference {
        name: "dpoly_degree_17_below",
        a: f64::from_bits(0x3fbd_180a_b612_cbde),
        b: f64::from_bits(0x3f59_f4df_9178_2eb9),
        root_hi: f64::from_bits(0x3fb2_07ab_4a58_326a),
        root_lo: f64::from_bits(0xbc33_2ad4_0240_3fc1),
        eta: f64::from_bits(0x3cb3_5196_acdb_316d),
    },
    Reference {
        name: "dpoly_degree_17_exact",
        a: f64::from_bits(0x3fbd_180a_b612_cbde),
        b: f64::from_bits(0x3f59_f4df_9178_2eba),
        root_hi: f64::from_bits(0x3fb2_07ab_4a58_326a),
        root_lo: f64::from_bits(0x3c31_b772_2f0d_68ae),
        eta: f64::from_bits(0x3cb3_5196_acdb_316d),
    },
    Reference {
        name: "dpoly_degree_17_above",
        a: f64::from_bits(0x3fbd_180a_b612_cbde),
        b: f64::from_bits(0x3f59_f4df_9178_2ebb),
        root_hi: f64::from_bits(0x3fb2_07ab_4a58_326a),
        root_lo: f64::from_bits(0x3c4b_4cdc_302d_888e),
        eta: f64::from_bits(0x3cb3_5196_acdb_316d),
    },
];

#[test]
fn reduced_degree_polynomials_retain_attainable_precision() {
    for r in REFERENCES {
        let s = Experimental::implied_total_volatility(-r.a, r.b).unwrap();
        assert!(s.is_finite() && s > 0.0, "{}: {s}", r.name);
        // The leading subtraction is exact by Sterbenz; the low component is
        // independent reference data, not a comparison with the old output.
        let rho = (((s - r.root_hi) - r.root_lo).abs() / r.root_hi) / r.eta;
        assert!(rho < 1.0, "{}: rho_J={rho}", r.name);
    }
}
