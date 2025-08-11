/// This example demonstrates how to use custom special functions with the `implied_vol` crate
/// by implementing the `SpecialFn` trait.
use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatility, PriceBlackScholes, SpecialFn};

#[link(name = "m")]
unsafe extern "C" {
    fn erf(x: f64) -> f64;
    fn erfc(x: f64) -> f64;
}

struct MySpecialFn;

// Implement the `SpecialFn` trait for `MySpecialFn`
// For more details on the `SpecialFn` trait, see the `implied_vol` crate documentation.
impl SpecialFn for MySpecialFn {
    fn erf(x: f64) -> f64 {
        // Call the C library function `erf` directly
        unsafe { erf(x) }
    }

    fn erfc(x: f64) -> f64 {
        // Call the C library function `erfc` directly
        unsafe { erfc(x) }
    }
}

fn main() {
    let forward = 100.0;
    let strike = 100.0;
    let maturity = 1.0;
    let volatility = 0.2;
    let is_call = true;

    let price_builder = PriceBlackScholes::builder()
        .forward(forward)
        .strike(strike)
        .expiry(maturity)
        .is_call(is_call)
        .vol(volatility)
        .build()
        .unwrap();

    let implied_vol_builder = ImpliedBlackVolatility::builder()
        .forward(forward)
        .strike(strike)
        .expiry(maturity)
        .is_call(is_call);

    // Use the custom special function implementation
    let price = price_builder.calculate::<MySpecialFn>().unwrap();
    let implied_vol = implied_vol_builder
        .option_price(price)
        .build()
        .unwrap()
        .calculate::<MySpecialFn>()
        .unwrap();
    println!("Price: {price}");
    println!("Implied Volatility: {implied_vol}");

    // Using the default special function implementation
    let price = price_builder.calculate::<DefaultSpecialFn>().unwrap();
    let implied_vol_builder = ImpliedBlackVolatility::builder()
        .forward(forward)
        .strike(strike)
        .expiry(maturity)
        .is_call(is_call);
    let implied_vol = implied_vol_builder
        .option_price(price)
        .build()
        .unwrap()
        .calculate::<DefaultSpecialFn>()
        .unwrap();
    println!("Price: {price}");
    println!("Implied Volatility: {implied_vol}");
}
