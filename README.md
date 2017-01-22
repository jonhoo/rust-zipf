# rust-zipf

[![crates.io page](https://img.shields.io/crates/v/zipf.svg)](https://crates.io/crates/zipf)
[![Downloads](https://img.shields.io/crates/d/zipf.png)](https://crates.io/crates/zipf)
[![Documentation](https://docs.rs/zipf/badge.svg)](https://docs.rs/zipf/)
[![Build Status](https://travis-ci.org/jonhoo/rust-zipf.svg?branch=master)](https://travis-ci.org/jonhoo/rust-zipf)

Rust implementation of a fast, discrete, bounded,
[Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law) random
number generator

For a random variable `X` whose values are distributed according to this distribution, the
probability mass function is given by

```ignore
P(X = k) = H(N,s) * 1 / k^s for k = 1,2,...,N
```

`H(N,s)` is the normalizing constant which corresponds to the generalized harmonic number
of order `N` of `s`.


This implementation is effectively a direct port of Apache Common's
[RejectionInversionZipfSampler](https://github.com/apache/commons-rng/blob/6a1b0c16090912e8fc5de2c1fb5bd8490ac14699/commons-rng-sampling/src/main/java/org/apache/commons/rng/sampling/distribution/RejectionInversionZipfSampler.java),
written in Java. It is based on the method described by Wolfgang Hörmann and Gerhard Derflinger
in [*Rejection-inversion to generate variates from monotone discrete
distributions*](https://dl.acm.org/citation.cfm?id=235029) from *ACM Transactions on Modeling
and Computer Simulation (TOMACS) 6.3 (1996)*.
