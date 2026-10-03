//! Exact-input Black roots at AS1 seed, geometry, quadrature, and scaling seams.
//! References use direct high-precision Black evaluation at two independent
//! precisions per input (at least 120 and 200 decimal digits). Inputs are exact
//! binary64 integer ratios; microscopic inputs use higher precision. The root
//! low part retains the difference from its nearest binary64 value.
//! These finite checks enforce the existing `rho_J < 1` core contract.
#![cfg(feature = "experimental")]

use implied_vol::ImpliedBlackVolatilityNormalised;
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
        name: "old_worst",
        a: f64::from_bits(0x3eb0_c6f7_a0b5_ed8d),
        b: f64::from_bits(0x3b73_0172_106d_797f),
        root_hi: f64::from_bits(0x3e81_aa07_2dd3_b703),
        root_lo: f64::from_bits(0x3b15_6a3d_d0ca_aaee),
        eta: f64::from_bits(0x3cb0_438b_8202_85e9),
    },
    Reference {
        name: "new_worst_and_largest_worsening",
        a: f64::from_bits(0x4020_0000_0000_0000),
        b: f64::from_bits(0x3ce2_200c_ebb2_a732),
        root_hi: f64::from_bits(0x3ff0_e229_8dcc_bf2f),
        root_lo: f64::from_bits(0x3c9e_1a88_266d_93b6),
        eta: f64::from_bits(0x3cb0_4420_c162_8cef),
    },
    Reference {
        name: "old_accepted_now_declined",
        a: f64::from_bits(0x4022_3333_3333_3333),
        b: f64::from_bits(0x3ce1_4bdd_4759_8fca),
        root_hi: f64::from_bits(0x3ff3_295c_1719_9e15),
        root_lo: f64::from_bits(0xbc9f_f200_368c_fbdd),
        eta: f64::from_bits(0x3cb0_43eb_2a0e_d0f9),
    },
    Reference {
        name: "rule12_transition",
        a: f64::from_bits(0x2b2b_ff2e_e48e_0531),
        b: f64::from_bits(0x23ec_489e_ba35_552d),
        root_hi: f64::from_bits(0x2af2_aa18_01d2_5499),
        root_lo: f64::from_bits(0x279f_d5d6_2954_cc30),
        eta: f64::from_bits(0x3cb0_1bdf_0476_887e),
    },
    Reference {
        name: "rule20_transition",
        a: f64::from_bits(0x3ff0_0000_0000_0000),
        b: f64::from_bits(0x2cf5_93e3_e893_012d),
        root_hi: f64::from_bits(0x3fa9_916c_1bce_9f77),
        root_lo: f64::from_bits(0x3c4f_aa3d_ebcb_4e7a),
        eta: f64::from_bits(0x3cb0_0a23_9532_35b9),
    },
    Reference {
        name: "scaled_rule1",
        a: f64::from_bits(0x2b2b_ff2e_e48e_052f),
        b: f64::from_bits(0x23ec_3757_6364_3d57),
        root_hi: f64::from_bits(0x2af2_aa04_1f2e_2960),
        root_lo: f64::from_bits(0x2780_d9da_afbb_8678),
        eta: f64::from_bits(0x3cb0_1bde_ca41_7b27),
    },
    Reference {
        name: "scaled_rule2",
        a: f64::from_bits(0x2b2b_ff2e_e48e_052f),
        b: f64::from_bits(0x183f_3ab5_d79a_ee06),
        root_hi: f64::from_bits(0x2ae6_65bc_5ae1_e7f0),
        root_lo: f64::from_bits(0xa78e_a0c3_9738_5b25),
        eta: f64::from_bits(0x3cb0_0a2a_02d3_b743),
    },
    Reference {
        name: "scaled_rule3",
        a: f64::from_bits(0x1668_7e92_154e_f7ac),
        b: f64::from_bits(0x1320_d477_9e7b_c1ab),
        root_hi: f64::from_bits(0x1639_9448_8a49_a796),
        root_lo: f64::from_bits(0x12c4_2fa1_736d_1500),
        eta: f64::from_bits(0x3cb0_427f_c06d_7f8c),
    },
    Reference {
        name: "dispatch_below_tiny_a_below",
        a: f64::from_bits(0x2b2b_ff2e_e48e_052f),
        b: f64::from_bits(0x27ef_b72a_c4c5_279d),
        root_hi: f64::from_bits(0x2afd_7a16_abce_6684),
        root_lo: f64::from_bits(0xa78b_cfc0_4dd1_d6b9),
        eta: f64::from_bits(0x3cb0_438b_8214_5723),
    },
    Reference {
        name: "dispatch_at_tiny_a_below",
        a: f64::from_bits(0x2b2b_ff2e_e48e_052f),
        b: f64::from_bits(0x27ef_b72a_c4c5_279e),
        root_hi: f64::from_bits(0x2afd_7a16_abce_6684),
        root_lo: f64::from_bits(0xa789_d988_3020_9665),
        eta: f64::from_bits(0x3cb0_438b_8214_5723),
    },
    Reference {
        name: "dispatch_above_tiny_a_below",
        a: f64::from_bits(0x2b2b_ff2e_e48e_052f),
        b: f64::from_bits(0x27ef_b72a_c4c5_279f),
        root_hi: f64::from_bits(0x2afd_7a16_abce_6684),
        root_lo: f64::from_bits(0xa787_e350_126f_5611),
        eta: f64::from_bits(0x3cb0_438b_8214_5723),
    },
    Reference {
        name: "dispatch_below_a_cutoff_equal",
        a: f64::from_bits(0x2b2b_ff2e_e48e_0530),
        b: f64::from_bits(0x27ef_b72a_c4c5_279e),
        root_hi: f64::from_bits(0x2afd_7a16_abce_6685),
        root_lo: f64::from_bits(0xa785_4e1f_51d1_d57d),
        eta: f64::from_bits(0x3cb0_438b_8214_5723),
    },
    Reference {
        name: "dispatch_below_a_cutoff_above",
        a: f64::from_bits(0x2b2b_ff2e_e48e_0531),
        b: f64::from_bits(0x27ef_b72a_c4c5_279f),
        root_hi: f64::from_bits(0x2afd_7a16_abce_6686),
        root_lo: f64::from_bits(0xa77d_98fc_aba3_a882),
        eta: f64::from_bits(0x3cb0_438b_8214_5723),
    },
    Reference {
        name: "deep_normal_price",
        a: f64::from_bits(0x4008_0000_0000_0000),
        b: f64::from_bits(0x0071_dde7_9e25_ac5c),
        root_hi: f64::from_bits(0x3fb4_9f52_4323_0edd),
        root_lo: f64::from_bits(0xbc5a_be25_4462_80a6),
        eta: f64::from_bits(0x3cb0_02f2_6d4b_4a1b),
    },
    Reference {
        name: "tiny_normal_price",
        a: f64::from_bits(0x03b8_f2b0_61ae_a072),
        b: f64::from_bits(0x0071_2440_5c42_d4fa),
        root_hi: f64::from_bits(0x038a_0d8b_5d65_bbdd),
        root_lo: f64::from_bits(0x0011_4e0c_82cb_695d),
        eta: f64::from_bits(0x3cb0_427f_c06d_7f8c),
    },
    Reference {
        name: "high_a_guard_fallback",
        a: f64::from_bits(0x4023_0000_0000_0000),
        b: f64::from_bits(0x3cdd_5037_0ebb_c4b7),
        root_hi: f64::from_bits(0x3ff3_f0a8_9f73_5090),
        root_lo: f64::from_bits(0x3c98_a72c_f67d_5f0a),
        eta: f64::from_bits(0x3cb0_4389_0bda_15d4),
    },
];

#[test]
fn independent_roots_cover_as1_seed_and_routing_changes() {
    for r in REFERENCES {
        let direct = Experimental::implied_total_volatility(-r.a, r.b).unwrap();
        assert!(direct.is_finite() && direct > 0.0, "{}: {direct}", r.name);
        // Sterbenz makes the leading subtraction exact. The independent low
        // part keeps this a comparison with the mathematical root.
        let rho = (((direct - r.root_hi) - r.root_lo).abs() / r.root_hi) / r.eta;
        assert!(rho < 1.0, "{}: rho_J={rho}", r.name);
        for x in [-r.a, r.a] {
            let public = ImpliedBlackVolatilityNormalised::builder()
                .log_moneyness(x)
                .normalised_price(r.b)
                .build()
                .unwrap()
                .calculate_with::<Experimental>()
                .unwrap();
            assert_eq!(public.to_bits(), direct.to_bits(), "{}: x={x}", r.name);
        }
    }
}
