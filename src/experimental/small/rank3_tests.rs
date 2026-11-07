//! Bit-for-bit comparison with the preceding SIMD factor layout.
use super::*;

// Frozen evaluator from before rank-three x/y repacking. Each live
// factor must retain its own fused Horner chain and scalar sum order.
fn original_layout<
    const X_DEGREE: usize,
    const Y_DEGREE: usize,
    const LANES: usize,
    const LENGTH: usize,
>(
    c: &[f64; LENGTH],
    rank: usize,
    x: f64,
    y: f64,
) -> f64 {
    let half = LANES / 2;
    let mut values = [0.0; 8];
    for k in (0..LANES).step_by(2) {
        let z = F64x2::splat(if k < half { x } else { y });
        let n = if k < half { X_DEGREE } else { Y_DEGREE };
        let mut v = F64x2::new(c[n * LANES + k], c[n * LANES + k + 1]);
        for j in (0..n).rev() {
            v = v.mul_add(z, F64x2::new(c[j * LANES + k], c[j * LANES + k + 1]));
        }
        values[k..k + 2].copy_from_slice(&v.to_array());
    }
    let mut result = values[0] * values[half];
    for k in 1..rank {
        result = fma(values[k], values[half + k], result);
    }
    result
}

fn check_table<
    const X_DEGREE: usize,
    const Y_DEGREE: usize,
    const LANES: usize,
    const LENGTH: usize,
>(
    c: &[f64; LENGTH],
    rank: usize,
) {
    assert_eq!(rank, 3);
    let inputs = [
        -1.0,
        -0.999_999_999_999_999_9,
        -0.75,
        -0.5,
        -0.1,
        -f64::EPSILON,
        -0.0,
        0.0,
        f64::EPSILON,
        0.1,
        0.5,
        0.75,
        0.999_999_999_999_999_9,
        1.0,
    ];
    for &x in &inputs {
        for &y in &inputs {
            let got = rank_seed_shape::<X_DEGREE, Y_DEGREE, LANES, LENGTH>(c, rank, x, y);
            let expected = original_layout::<X_DEGREE, Y_DEGREE, LANES, LENGTH>(c, rank, x, y);
            assert_eq!(
                got.to_bits(),
                expected.to_bits(),
                "degrees={X_DEGREE}/{Y_DEGREE}, x={x}, y={y}"
            );
        }
    }
}

#[test]
fn all_rank_three_seed_tables_preserve_original_bits() {
    for (index, &(c, rank)) in small_lr.iter().enumerate() {
        if rank != 3 {
            continue;
        }
        match SMALL_RANK_DEGREES[index] {
            [6, 5] => check_table::<6, 5, 8, 56>(c.first_chunk::<56>().unwrap(), rank),
            [7, 5] => check_table::<7, 5, 8, 64>(c.first_chunk::<64>().unwrap(), rank),
            degrees => panic!("unexpected small degrees: {degrees:?}"),
        }
    }
    for (index, &(c, rank)) in core_lr.iter().enumerate() {
        if rank != 3 {
            continue;
        }
        match CORE_RANK_DEGREES[index] {
            [6, 5] => check_table::<6, 5, 8, 56>(c.first_chunk::<56>().unwrap(), rank),
            [6, 6] => check_table::<6, 6, 8, 56>(c.first_chunk::<56>().unwrap(), rank),
            [7, 5] => check_table::<7, 5, 8, 64>(c.first_chunk::<64>().unwrap(), rank),
            [7, 6] => check_table::<7, 6, 8, 64>(c.first_chunk::<64>().unwrap(), rank),
            degrees => panic!("unexpected core degrees: {degrees:?}"),
        }
    }
    for &(c, rank) in &factors {
        if rank == 3 {
            check_table::<8, 8, 8, 72>(c.first_chunk::<72>().unwrap(), rank);
        }
    }
}
