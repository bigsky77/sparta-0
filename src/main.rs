#[macro_use]
extern crate lazy_static;
extern crate nalgebra as na;
extern crate rand;
extern crate num_bigint;
extern crate num_traits;

mod polynomial;
mod prover;
mod verifier;

use prover::prove;
use verifier::verify;

use polynomial::{FiniteField, PolynomialRing};
use na::DVector;

fn main() {

    println!();
    println!("\x1b[32mSparta(0) Starting Execution\x1b[0m");

    // Finite field prime
    let prime = 7;

    // Variable names as a single string
    let variable_names = "x1, x2, y1, y2, y3, x11, x22, y11, y22, y33";

    // Create PolynomialRing
    let ring = PolynomialRing::new(prime, 10, variable_names);

    // Generate variables
    let _vars = ring.gens();

    let x = vec![0, 1, 1, 2, 3, 6, 6];
    let w: Vec<i32> = Vec::new();

    // set values to finite field
    let combined: Vec<i32> = x.iter().chain(&w).cloned().collect();
    let combined_ff: Vec<FiniteField> = combined.iter().map(|&v| FiniteField::new(v, prime)).collect();

    let mut extended_combined_ff = combined_ff.clone();
    extended_combined_ff.push(FiniteField::new(1, prime));

    let combined_values: Vec<i32> = extended_combined_ff.iter().map(|ff| ff.get_value()).collect();
    let z1 = DVector::from_vec(combined_values);
    println!("      Z1: {:?}", z1);

    // Run the prover
    let proof = prove(z1);

    // Run the verifier
    verify(proof);

    println!();
    println!("\x1b[32mSparta(0) Execution Complete\x1b[0m");
}
