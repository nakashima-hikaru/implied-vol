pub(crate) const SQRT_PI_OVER_TWO: f64 = f64::from_bits(0x3ff40d931ff62705);
pub(crate) const SQRT_TWO_PI: f64 = f64::from_bits(0x40040d931ff62705);
pub(crate) const FRAC_SQRT_TWO_PI: f64 = f64::from_bits(0x3fd9884533d43651);
pub(crate) const SQRT_THREE: f64 = f64::from_bits(0x3ffbb67ae8584caa);
pub(crate) const SQRT_TWO_OVER_PI: f64 = f64::from_bits(0x3fe9884533d43651);
pub(crate) const ONE_OVER_SQRT_THREE: f64 = f64::from_bits(0x3fe279a74590331d);
pub(crate) const TWO_PI_OVER_SQRT_TWENTY_SEVEN: f64 = f64::from_bits(0x3ff358e1a79ed7e1);
pub(crate) const SQRT_THREE_OVER_THIRD_ROOT_TWO_PI: f64 = f64::from_bits(0x3fee095e112353f7);
pub(crate) const FOURTH_ROOT_DBL_EPSILON: f64 = 0.0001220703125;
pub(crate) const SIXTEENTH_ROOT_DBL_EPSILON: f64 = 0.10511205190671433;
pub(crate) const SQRT_MIN_POSITIVE: f64 = 1.4916681462400413e-154;
pub(crate) const SQRT_DBL_MAX: f64 = 1.3407807929942596e154;
// pub(crate) const SQRT_DBL_MIN: f64 = 1.3407807929942596e154;
// Set this to 0 if you want positive results for (positive) denormalised inputs, else to DBL_MIN.
// Note that you cannot achieve full machine accuracy from denormalised inputs!
pub(crate) const DENORMALISATION_CUTOFF: f64 = 0.0;
pub(crate) const ONE_OVER_SQRT_TWO_PI: f64 = f64::from_bits(0x3fd9884533d43651);
pub(crate) const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC: f64 = f64::NEG_INFINITY;
pub(crate) const VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM: f64 = f64::INFINITY;
pub(crate) const HALF_OF_LN_TWO_PI: f64 = f64::from_bits(0x3fed67f1c864beb4);
