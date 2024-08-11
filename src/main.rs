#[macro_use]
extern crate lazy_static;
extern crate nalgebra as na;
extern crate rand;
extern crate num_bigint;
extern crate num_traits;
extern crate ronkathon;

use std::collections::BTreeMap;
use std::ops::{Mul, Sub};

use ronkathon::polynomial::multivariate_polynomial::{MultivariateTerm, MultivariateVariable};
use ronkathon::{polynomial::multivariate_polynomial::MultivariatePolynomial, algebra::field::prime::PlutoBaseField, algebra::field::prime::PlutoPrime, algebra::field::prime::PrimeField};
use ronkathon::sumcheck::{non_interactive_sumcheck_prove, non_interactive_sumcheck_verify};
use ronkathon::random::Random;

use rand::{Rng, random};
use na::{DVector, DMatrix};

lazy_static! {
    pub static ref M1: DMatrix<PlutoBaseField> = DMatrix::from_row_slice(4, 8, &[
        PlutoBaseField::from(1), PlutoBaseField::from(1), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(1), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(1), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0),
   ]);

    pub static ref M2: DMatrix<PlutoBaseField> = DMatrix::from_row_slice(4, 8, &[
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(1),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(1), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(1), PlutoBaseField::from(0), PlutoBaseField::from(0),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0),
   ]);

    pub static ref M3: DMatrix<PlutoBaseField> = DMatrix::from_row_slice(4, 8, &[
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(1), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(1), PlutoBaseField::from(0), PlutoBaseField::from(0),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(1), PlutoBaseField::from(0),
        PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0), PlutoBaseField::from(0),
   ]);
}

fn main() {

    println!();
    println!("\x1b[32mSparta(0): Starting Execution\x1b[0m");

    let poly = MultivariatePolynomial::<PlutoBaseField>::from_terms(vec![
            MultivariateTerm::new(
            vec![MultivariateVariable { index: 1, exponent: 1 }],
            PlutoBaseField::new(1)
            ),
            MultivariateTerm::new(
            vec![MultivariateVariable { index: 2, exponent: 1 }],
            PlutoBaseField::new(1)
            ),
            MultivariateTerm::new(
            vec![MultivariateVariable { index: 3, exponent: 1 }],
            PlutoBaseField::new(2)
            ),
            MultivariateTerm::new(
            vec![MultivariateVariable { index: 4, exponent: 1 }],
            PlutoBaseField::new(3)
            ),
            MultivariateTerm::new(
            vec![MultivariateVariable { index: 5, exponent: 1 }],
            PlutoBaseField::new(6)
            ),
            MultivariateTerm::new(
            vec![MultivariateVariable { index: 6, exponent: 1 }],
            PlutoBaseField::new(6)
            ),
            MultivariateTerm::new(
            vec![MultivariateVariable { index: 0, exponent: 1 }],
            PlutoBaseField::new(1)
            )
        ]);

    println!("Z1 Polynomial: {:?}", poly.to_string());

    let mz1 = &poly.clone().mul_matrix(&M1);
    let mz2 = &poly.clone().mul_matrix(&M2);
    let mz3 = &poly.clone().mul_matrix(&M3);

    let elementwise_product = &mz1.clone().mul(mz2.clone());
    let z_poly = elementwise_product.clone().sub(mz3.clone());

    //// Run the prover
    println!();
    println!("\x1b[32mSparta(0): Running Prover\x1b[0m");
    let proof = non_interactive_sumcheck_prove::<PlutoBaseField>(&z_poly);

    //// Run the verifier
    println!();
    println!("\x1b[32mSparta(0): Running Verifier\x1b[0m");
    non_interactive_sumcheck_verify(&proof, &z_poly);

    println!();
    println!("\x1b[32mSparta(0): Execution Complete\x1b[0m");
}
