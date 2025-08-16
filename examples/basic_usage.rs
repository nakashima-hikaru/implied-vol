use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatility};

fn main() {
    let implied_vol_builder = ImpliedBlackVolatility::builder()
        .option_price(10.0)
        .forward(100.0)
        .strike(100.0)
        .expiry(1.0)
        .is_call(true)
        .build();

    // Check if the inputs of builder are valid:
    // - Option price must be non-negative and finite (e.g., not NaN or Inf)
    // - Forward must be non-negative and finite
    // - Strike must be non-negative and finite
    // - Expiry must be non-negative but may be infinite
    assert!(implied_vol_builder.is_some());
    let implied_vol = implied_vol_builder.unwrap().calculate::<DefaultSpecialFn>();
    // If None, then implied volatility does not exist under the given inputs:
    assert!(implied_vol.is_some());
    println!("{}", implied_vol.unwrap());
}
