#![feature(test)]
extern crate test;
extern crate zipf;
extern crate rand;
extern crate randomkit;

use test::Bencher;

#[bench]
fn bench_us(b: &mut Bencher) {
    use rand;
    use rand::distributions::Sample;
    use zipf::ZipfDistribution;
    let mut rng = rand::thread_rng();
    let mut us = ZipfDistribution::new(1000000, 1.07).unwrap();
    b.iter(|| us.sample(&mut rng));
}

#[bench]
fn bench_randomkit(b: &mut Bencher) {
    use randomkit::Sample;
    let mut rng = randomkit::Rng::new().unwrap();
    let oracle = randomkit::dist::Zipf::new(1.07).unwrap();
    b.iter(|| oracle.sample(&mut rng));
}

#[bench]
fn bench_threadrng(b: &mut Bencher) {
    use rand::{self, Rng};
    let mut rng = rand::thread_rng();
    b.iter(|| rng.gen::<u64>());
}
