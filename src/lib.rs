//! A fast generator of discrete, bounded
//! [Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law) random numbers.
//!
//! For a random variable `X` whose values are distributed according to this distribution, the
//! probability mass function is given by
//!
//! ```ignore
//! P(X = k) = H(N,s) * 1 / k^s for k = 1,2,...,N
//! ```
//!
//! `H(N,s)` is the normalizing constant which corresponds to the generalized harmonic number
//! of order `N` of `s`.
//!
//! # Example
//!
//! ```
//! use rand::distributions::Distribution;
//!
//! let mut rng = rand::thread_rng();
//! let mut zipf = zipf::ZipfDistribution::new(1000, 1.03).unwrap();
//! let sample = zipf.sample(&mut rng);
//! ```
//!
//! This implementation is effectively a direct port of Apache Common's
//! [RejectionInversionZipfSampler][ref], written in Java. It is based on the method described by
//! Wolfgang Hörmann and Gerhard Derflinger in [*Rejection-inversion to generate variates from
//! monotone discrete distributions*](https://dl.acm.org/citation.cfm?id=235029) from *ACM
//! Transactions on Modeling and Computer Simulation (TOMACS) 6.3 (1996)*.
//!
//! [ref]: https://github.com/apache/commons-rng/blob/6a1b0c16090912e8fc5de2c1fb5bd8490ac14699/commons-rng-sampling/src/main/java/org/apache/commons/rng/sampling/distribution/RejectionInversionZipfSampler.java

#![warn(rust_2018_idioms)]

use rand::Rng;

/// Random number generator that generates Zipf-distributed random numbers using rejection
/// inversion.
#[derive(Clone, Copy)]
pub struct ZipfDistribution {
    /// Number of elements
    num_elements: f64,
    /// Exponent parameter of the distribution
    exponent: f64,
    /// `hIntegral(1.5) - 1}`
    h_integral_x1: f64,
    /// `hIntegral(num_elements + 0.5)}`
    h_integral_num_elements: f64,
    /// `2 - hIntegralInverse(hIntegral(2.5) - h(2)}`
    s: f64,
}

impl ZipfDistribution {
    /// Creates a new [Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law)
    /// random number generator.
    ///
    /// Note that both the number of elements and the exponent must be greater than 0.
    pub fn new(num_elements: usize, exponent: f64) -> Result<Self, ()> {
        if num_elements == 0 {
            return Err(());
        }
        if exponent <= 0f64 {
            return Err(());
        }

        let z = ZipfDistribution {
            num_elements: num_elements as f64,
            exponent,
            h_integral_x1: ZipfDistribution::h_integral(1.5, exponent) - 1f64,
            h_integral_num_elements: ZipfDistribution::h_integral(
                num_elements as f64 + 0.5,
                exponent,
            ),
            s: 2f64
                - ZipfDistribution::h_integral_inv(
                    ZipfDistribution::h_integral(2.5, exponent)
                        - ZipfDistribution::h(2f64, exponent),
                    exponent,
                ),
        };

        // populate cache

        Ok(z)
    }
}

impl ZipfDistribution {
    fn next<R: Rng + ?Sized>(&self, rng: &mut R) -> usize {
        // The paper describes an algorithm for exponents larger than 1 (Algorithm ZRI).
        //
        // The original method uses
        //   H(x) = (v + x)^(1 - q) / (1 - q)
        // as the integral of the hat function.
        //
        // This function is undefined for q = 1, which is the reason for the limitation of the
        // exponent.
        //
        // If instead the integral function
        //   H(x) = ((v + x)^(1 - q) - 1) / (1 - q)
        // is used, for which a meaningful limit exists for q = 1, the method works for all
        // positive exponents.
        //
        // The following implementation uses v = 0 and generates integral number in the range [1,
        // num_elements]. This is different to the original method where v is defined to
        // be positive and numbers are taken from [0, i_max]. This explains why the implementation
        // looks slightly different.

        let hnum = self.h_integral_num_elements;

        loop {
            use std::cmp;
            let u: f64 = hnum + rng.gen::<f64>() * (self.h_integral_x1 - hnum);
            // u is uniformly distributed in (h_integral_x1, h_integral_num_elements]

            let x: f64 = ZipfDistribution::h_integral_inv(u, self.exponent);

            // Limit k to the range [1, num_elements] if it would be outside
            // due to numerical inaccuracies.
            let k64 = x.max(1.0).min(self.num_elements);
            // float -> integer rounds towards zero
            let k = cmp::max(1, k64 as usize);

            // Here, the distribution of k is given by:
            //
            //   P(k = 1) = C * (hIntegral(1.5) - h_integral_x1) = C
            //   P(k = m) = C * (hIntegral(m + 1/2) - hIntegral(m - 1/2)) for m >= 2
            //
            // where C = 1 / (h_integral_num_elements - h_integral_x1)
            if k64 - x <= self.s
                || u >= ZipfDistribution::h_integral(k64 + 0.5, self.exponent)
                    - ZipfDistribution::h(k64, self.exponent)
            {
                // Case k = 1:
                //
                //   The right inequality is always true, because replacing k by 1 gives
                //   u >= hIntegral(1.5) - h(1) = h_integral_x1 and u is taken from
                //   (h_integral_x1, h_integral_num_elements].
                //
                //   Therefore, the acceptance rate for k = 1 is P(accepted | k = 1) = 1
                //   and the probability that 1 is returned as random value is
                //   P(k = 1 and accepted) = P(accepted | k = 1) * P(k = 1) = C = C / 1^exponent
                //
                // Case k >= 2:
                //
                //   The left inequality (k - x <= s) is just a short cut
                //   to avoid the more expensive evaluation of the right inequality
                //   (u >= hIntegral(k + 0.5) - h(k)) in many cases.
                //
                //   If the left inequality is true, the right inequality is also true:
                //     Theorem 2 in the paper is valid for all positive exponents, because
                //     the requirements h'(x) = -exponent/x^(exponent + 1) < 0 and
                //     (-1/hInverse'(x))'' = (1+1/exponent) * x^(1/exponent-1) >= 0
                //     are both fulfilled.
                //     Therefore, f(x) = x - hIntegralInverse(hIntegral(x + 0.5) - h(x))
                //     is a non-decreasing function. If k - x <= s holds,
                //     k - x <= s + f(k) - f(2) is obviously also true which is equivalent to
                //     -x <= -hIntegralInverse(hIntegral(k + 0.5) - h(k)),
                //     -hIntegralInverse(u) <= -hIntegralInverse(hIntegral(k + 0.5) - h(k)),
                //     and finally u >= hIntegral(k + 0.5) - h(k).
                //
                //   Hence, the right inequality determines the acceptance rate:
                //   P(accepted | k = m) = h(m) / (hIntegrated(m+1/2) - hIntegrated(m-1/2))
                //   The probability that m is returned is given by
                //   P(k = m and accepted) = P(accepted | k = m) * P(k = m)
                //                         = C * h(m) = C / m^exponent.
                //
                // In both cases the probabilities are proportional to the probability mass
                // function of the Zipf distribution.

                return k;
            }
        }
    }
}

impl rand::distributions::Distribution<usize> for ZipfDistribution {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> usize {
        self.next(rng)
    }
}

use std::fmt;
impl fmt::Debug for ZipfDistribution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        f.debug_struct("ZipfDistribution")
            .field("e", &self.exponent)
            .field("n", &self.num_elements)
            .finish()
    }
}

impl ZipfDistribution {
    /// Computes `H(x)`, defined as
    ///
    ///  - `(x^(1 - exponent) - 1) / (1 - exponent)`, if `exponent != 1`
    ///  - `log(x)`, if `exponent == 1`
    ///
    /// `H(x)` is an integral function of `h(x)`, the derivative of `H(x)` is `h(x)`.
    fn h_integral(x: f64, exponent: f64) -> f64 {
        let log_x = x.ln();
        helper2((1f64 - exponent) * log_x) * log_x
    }

    /// Computes `h(x) = 1 / x^exponent`
    fn h(x: f64, exponent: f64) -> f64 {
        (-exponent * x.ln()).exp()
    }

    /// The inverse function of `H(x)`.
    /// Returns the `y` for which `H(y) = x`.
    fn h_integral_inv(x: f64, exponent: f64) -> f64 {
        let mut t: f64 = x * (1f64 - exponent);
        if t < -1f64 {
            // Limit value to the range [-1, +inf).
            // t could be smaller than -1 in some rare cases due to numerical errors.
            t = -1f64;
        }
        (helper1(t) * x).exp()
    }
}

/// Helper function that calculates `log(1 + x) / x`.
/// A Taylor series expansion is used, if x is close to 0.
fn helper1(x: f64) -> f64 {
    if x.abs() > 1e-8 {
        x.ln_1p() / x
    } else {
        1f64 - x * (0.5 - x * (1.0 / 3.0 - 0.25 * x))
    }
}

/// Helper function to calculate `(exp(x) - 1) / x`.
/// A Taylor series expansion is used, if x is close to 0.
fn helper2(x: f64) -> f64 {
    if x.abs() > 1e-8 {
        x.exp_m1() / x
    } else {
        1f64 + x * 0.5 * (1f64 + x * 1.0 / 3.0 * (1f64 + 0.25 * x))
    }
}

#[cfg(test)]
#[cfg_attr(tarpaulin, skip)]
mod test {
    use super::ZipfDistribution;
    use rand::distributions::Distribution;

    #[inline]
    fn test(alpha: f64) {
        const N: usize = 100;

        // as the alpha increases, we need more samples, since the frequency in the tail grows so low
        let samples = (2.0f64.powf(alpha) * 5000000.0) as usize;

        let mut rng = rand::thread_rng();
        let zipf = ZipfDistribution::new(N, alpha).unwrap();

        let harmonic: f64 = (1..=N).map(|n| 1.0 / (n as f64).powf(alpha)).sum();

        // sample a bunch
        let mut buckets = vec![0; N];
        for _ in 0..samples {
            let sample = zipf.sample(&mut rng);
            buckets[sample - 1] += 1;
        }

        // for each bucket, see that frequency roughly matches what we expect
        // note that the ratios here are ratios _of fractions_, so 0.5 does not mean we're off by
        // 50%, it means we're off by 50% _of the frequency_. in other words, if the frequency was
        // supposed to be 0.10, and it was actually 0.05, that's a 50% difference.
        for i in 0..N {
            let freq = buckets[i] as f64 / samples as f64;
            let expected = (1.0 / ((i + 1) as f64).powf(alpha)) / harmonic;
            // println!("{} {} {}", i + 1, freq, expected);

            let off_by = (expected - freq).abs();
            assert!(off_by < 0.1); // never be off by more than 10% in absolute terms

            let good = if i == 0 {
                // the first point tends to overshoot by a bit
                // that overshoot is spread out across the remaining point
                // which all tend to undershoot slightly
                freq > expected && off_by < 0.3 * expected
            } else if i == N - 1 {
                // last one undershoots by the remainder
                // so it tends to be off by a fair amount
                // though it's relative, since the frequency is also so low
                off_by < expected
            } else {
                // others should only marginally undershoot
                off_by < 0.5 * expected
            };
            if !good {
                panic!("got {}, expected {} for k = {}", freq, expected, i + 1);
            }
        }
    }

    #[test]
    fn one() {
        test(1.00);
    }

    #[test]
    // until https://github.com/xd009642/tarpaulin/issues/486 lands
    // #[cfg_attr(tarpaulin, ignore)]
    #[ignore]
    fn two() {
        test(2.00);
    }

    #[test]
    // until https://github.com/xd009642/tarpaulin/issues/486 lands
    // #[cfg_attr(tarpaulin, ignore)]
    #[ignore]
    fn three() {
        test(3.00);
    }

    #[test]
    fn float() {
        test(1.08);
    }

    #[test]
    fn debug() {
        eprintln!("{:?}", ZipfDistribution::new(100, 1.0).unwrap());
    }

    #[test]
    fn errs() {
        ZipfDistribution::new(0, 1.0).unwrap_err();
        ZipfDistribution::new(100, 0.0).unwrap_err();
        ZipfDistribution::new(100, -1.0).unwrap_err();
    }
}
