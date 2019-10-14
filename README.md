# rust-zipf

[![Crates.io](https://img.shields.io/crates/v/zipf.svg)](https://crates.io/crates/zipf)
[![Documentation](https://docs.rs/zipf/badge.svg)](https://docs.rs/zipf/)
[![Build Status](https://dev.azure.com/jonhoo/jonhoo/_apis/build/status/zipf?branchName=master)](https://dev.azure.com/jonhoo/jonhoo/_build/latest?definitionId=14&branchName=master)
[![Codecov](https://codecov.io/github/jonhoo/rust-zipf/coverage.svg?branch=master)](https://codecov.io/gh/jonhoo/rust-zipf)
[![Dependency status](https://deps.rs/repo/github/jonhoo/rust-zipf/status.svg)](https://deps.rs/repo/github/jonhoo/rust-zipf)

Rust implementation of a fast, discrete, bounded,
[Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law) random
number generator. Compared to the implementation provided by
[`randomkit`](https://github.com/stygstra/rust-randomkit) (which binds
to NumPy's fork of RandomKit), this crate is approximately twice as
fast:

```console
$ cargo +nightly bench
test tests::bench_randomkit ... bench:         339 ns/iter (+/- 18)
test tests::bench_us        ... bench:          68 ns/iter (+/- 1)
test tests::bench_threadrng ... bench:          11 ns/iter (+/- 0)
```

It is also both driven by, and provides, a [Rust random number
generator](https://doc.rust-lang.org/rand/rand/trait.Rng.html).

This implementation is effectively a direct port of Apache Common's
[RejectionInversionZipfSampler](https://github.com/apache/commons-rng/blob/6a1b0c16090912e8fc5de2c1fb5bd8490ac14699/commons-rng-sampling/src/main/java/org/apache/commons/rng/sampling/distribution/RejectionInversionZipfSampler.java),
written in Java. It is based on the method described by Wolfgang Hörmann and Gerhard Derflinger
in [*Rejection-inversion to generate variates from monotone discrete
distributions*](https://dl.acm.org/citation.cfm?id=235029) from *ACM Transactions on Modeling
and Computer Simulation (TOMACS) 6.3 (1996)*.
