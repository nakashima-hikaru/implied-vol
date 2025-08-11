/// This example demonstrates how to use custom special functions with the `implied_vol` crate
/// by implementing the `SpecialFn` trait.
use implied_vol::{DefaultSpecialFn, SpecialFn, black};

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

    // Use the custom special function implementation
    let price =
        black::price::<MySpecialFn>(forward, strike, volatility, maturity, is_call).unwrap();
    let implied_vol =
        black::implied_volatility::<MySpecialFn>(price, forward, strike, maturity, is_call)
            .unwrap();
    println!("Price: {price}");
    println!("Implied Volatility: {implied_vol}");

    // Using the default special function implementation
    let price =
        black::price::<DefaultSpecialFn>(forward, strike, volatility, maturity, is_call).unwrap();
    let implied_vol =
        black::implied_volatility::<DefaultSpecialFn>(price, forward, strike, maturity, is_call)
            .unwrap();
    println!("Price: {price}");
    println!("Implied Volatility: {implied_vol}");
}
