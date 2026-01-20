# Implied Volatility Accuracy Report

## Overview

This report documents the results of a Monte Carlo simulation designed to verify the accuracy of the normalized implied volatility calculation implemented in `lets_be_rational`. The goal was to check if the relative error of the recalculated volatility falls within machine epsilon (`f64::EPSILON`).

## Methodology

A test function `check_relative_error_of_implied_vol` was added to `src/lets_be_rational.rs`. The test procedure is as follows:

1.  Generate 1,000,000 random pairs of `theta_x` (log-moneyness) and `s` (total volatility).
    *   `theta_x` range: `[-5.0, -0.0001]`
    *   `s` range: `[0.001, 5.0]`
2.  Calculate the theoretical normalized option price `beta` using `bs_option_price::normalised_black`.
3.  Attempt to recalculate `s` (implied volatility) from `beta` and `theta_x` using the `lets_be_rational` function.
4.  Measure the relative error: `|s_original - s_recalc| / s_original`.
5.  Count occurrences where the relative error exceeds `f64::EPSILON`.
6.  Count occurrences where the calculation panics (specifically due to `s` becoming non-positive during iteration).

## Results

**Sample Size:** 1,000,000

| Metric | Value |
| :--- | :--- |
| **Max Relative Error** (excluding panics) | `2.50398e-1` (approx. 25%) |
| **Percentage Exceeding Machine Epsilon** | `10.94%` (109,440 samples) |
| **Panic Rate** | `0.14%` (1,395 samples) |
| **Machine Epsilon** | `2.22045e-16` |

## Analysis

The findings indicate that the current implementation **does not** consistently satisfy the requirement of relative error being less than machine epsilon.

1.  **High Max Error:** The maximum relative error observed is significantly large (~0.25), indicating that for certain input combinations (likely in edge case regions of deep OTM or high/low volatility), the algorithm converges to a suboptimal solution or diverges significantly before stopping.
2.  **Frequent Epsilon Violation:** Approximately 11% of the random samples resulted in errors greater than machine epsilon.
3.  **Stability Issues:** The implementation panicked in roughly 0.14% of cases due to the volatility parameter `s` becoming non-positive during the Householder iteration, which suggests numerical instability in specific regions.

## Conclusion

The `lets_be_rational` implementation, as tested, does not strictly meet the "relative error <= machine epsilon" criteria for all inputs within the tested domain.
