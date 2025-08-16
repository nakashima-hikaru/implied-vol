use crate::{SpecialFn, lets_be_rational};
use bon::Builder;

#[derive(Builder)]
#[builder(derive(Clone, Debug))]
#[builder(finish_fn(vis = "", name = build_internal))]
pub struct ImpliedBlackVolatility {
    forward: f64,
    strike: f64,
    expiry: f64,
    is_call: bool,
    option_price: f64,
}

impl<S: implied_black_volatility_builder::IsComplete> ImpliedBlackVolatilityBuilder<S> {
    pub fn build(self) -> Option<ImpliedBlackVolatility> {
        let implied_black_volatility = self.build_internal();

        if matches!(
            implied_black_volatility.forward.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) || implied_black_volatility.forward.is_infinite()
        {
            return None;
        }
        if matches!(
            implied_black_volatility.strike.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) || implied_black_volatility.strike.is_infinite()
        {
            return None;
        }
        if matches!(
            implied_black_volatility.expiry.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) {
            return None;
        }
        if matches!(
            implied_black_volatility.option_price.partial_cmp(&0.0),
            Some(std::cmp::Ordering::Less) | None
        ) || implied_black_volatility.option_price.is_infinite()
        {
            return None;
        }
        Some(implied_black_volatility)
    }

    pub fn build_unchecked(self) -> ImpliedBlackVolatility {
        self.build_internal()
    }
}

impl ImpliedBlackVolatility {
    #[must_use]
    #[inline(always)]
    pub fn calculate<SpFn: SpecialFn>(&self) -> Option<f64> {
        if self.is_call {
            lets_be_rational::implied_black_volatility_input_unchecked::<SpFn, true>(
                self.option_price,
                self.forward,
                self.strike,
                self.expiry,
            )
        } else {
            lets_be_rational::implied_black_volatility_input_unchecked::<SpFn, false>(
                self.option_price,
                self.forward,
                self.strike,
                self.expiry,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::DefaultSpecialFn;
    use crate::builder::implied_black_volatility::ImpliedBlackVolatility;

    #[test]
    fn strike_anomaly() {
        for k in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let price = 100.0;
            let f = 100.0;
            let t = 1.0;
            const Q: bool = true;
            assert!(
                ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(f)
                    .strike(k)
                    .expiry(t)
                    .is_call(Q)
                    .build()
                    .is_none()
            );
        }
    }

    #[test]
    fn forward_anomaly() {
        for f in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let price = 100.0;
            let k = 100.0;
            let t = 1.0;
            const Q: bool = true;
            assert!(
                ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(f)
                    .strike(k)
                    .expiry(t)
                    .is_call(Q)
                    .build()
                    .is_none()
            );
        }
    }

    #[test]
    fn price_anomaly() {
        for price in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let f = 100.0;
            let t = 1.0;
            let k = 100.0;
            const Q: bool = true;
            assert!(
                ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(f)
                    .strike(k)
                    .expiry(t)
                    .is_call(Q)
                    .build()
                    .is_none()
            );
        }
    }

    #[test]
    fn time_anomaly() {
        for t in [f64::NAN, f64::NEG_INFINITY] {
            let price = 10.0;
            let f = 100.0;
            let k = 100.0;
            const Q: bool = true;
            assert!(
                ImpliedBlackVolatility::builder()
                    .option_price(price)
                    .forward(f)
                    .strike(k)
                    .expiry(t)
                    .is_call(Q)
                    .build()
                    .is_none()
            );
        }
    }

    #[test]
    fn time_inf() {
        let price = 10.0;
        let f = 100.0;
        let k = 100.0;
        let t = f64::INFINITY;
        const Q: bool = true;

        let vol = ImpliedBlackVolatility::builder()
            .option_price(price)
            .forward(f)
            .strike(k)
            .expiry(t)
            .is_call(Q)
            .build()
            .unwrap()
            .calculate::<DefaultSpecialFn>()
            .unwrap();
        assert_eq!(vol, 0.0);
    }
}
