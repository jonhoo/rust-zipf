#![feature(test)]
extern crate rand;
extern crate test;
extern crate zipf;

use test::Bencher;

#[bench]
fn bench_us(b: &mut Bencher) {
    use rand;
    use rand::distributions::Distribution;
    use zipf::ZipfDistribution;
    let mut rng = rand::thread_rng();
    let us = ZipfDistribution::new(1000000, 1.07).unwrap();
    b.iter(|| us.sample(&mut rng));
}

#[bench]
fn bench_threadrng(b: &mut Bencher) {
    use rand::{self, Rng};
    let mut rng = rand::thread_rng();
    b.iter(|| rng.gen::<u64>());
}
