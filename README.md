# rust-zipf

[![crates.io page](https://img.shields.io/crates/v/zipf.svg)](https://crates.io/crates/zipf)
[![Downloads](https://img.shields.io/crates/d/zipf.png)](https://crates.io/crates/zipf)
[![Documentation](https://docs.rs/zipf/badge.svg)](https://docs.rs/zipf/)
[![Build Status](https://travis-ci.org/jonhoo/rust-zipf.svg?branch=master)](https://travis-ci.org/jonhoo/rust-zipf)

Rust implementation of a fast, discrete, bounded,
[Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law) random
number generator. Compared to the implementation provided by
[`randomkit`](https://github.com/stygstra/rust-randomkit) (which binds
to NumPy's fork of RandomKit), this crate is approximately twice as
fast:

```console
$ cargo +nightly bench
test tests::bench_randomkit ... bench:         344 ns/iter (+/- 15)
test tests::bench_us        ... bench:          66 ns/iter (+/- 2)
test tests::bench_threadrng ... bench:           9 ns/iter (+/- 0)
```

It is also both driven by, and provides, a [Rust random number
generator](https://doc.rust-lang.org/rand/rand/trait.Rng.html).

This implementation is effectively a direct port of Apache Common's
[RejectionInversionZipfSampler](https://github.com/apache/commons-rng/blob/6a1b0c16090912e8fc5de2c1fb5bd8490ac14699/commons-rng-sampling/src/main/java/org/apache/commons/rng/sampling/distribution/RejectionInversionZipfSampler.java),
written in Java. It is based on the method described by Wolfgang Hörmann and Gerhard Derflinger
in [*Rejection-inversion to generate variates from monotone discrete
distributions*](https://dl.acm.org/citation.cfm?id=235029) from *ACM Transactions on Modeling
and Computer Simulation (TOMACS) 6.3 (1996)*.
