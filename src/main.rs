#[macro_use]
extern crate lazy_static;
extern crate nalgebra as na;
extern crate rand;
extern crate num_bigint;
extern crate num_traits;
extern crate ronkathon;

use std::collections::BTreeMap;

use ronkathon::{polynomial::multivariate_polynomial::MultivariatePolynomial, algebra::field::prime::PlutoBaseField, algebra::field::prime::PlutoPrime, algebra::field::prime::PrimeField};
use ronkathon::polynomial::sumcheck::non_interactive_sumcheck_prove;
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

impl Random for PlutoBaseField {
  fn random<R: Rng + ?Sized>(rng: &mut R) -> Self {
    let value = rng.gen_range(0..PlutoPrime::Base as usize);
    PlutoBaseField::new(value)
  }
}

fn main() {

    println!();
    println!("\x1b[32mSparta(0) Starting Execution\x1b[0m");

    // Variable names as a single string
    let variable_names = "x1, x2, y1, y2, y3, x11, x22, y11, y22, y33";

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
    //let mz1 = poly.mul_matrix(&M1);
    //let mz2 = poly.mul_matrix(&M1);
    //let mz3 = poly.mul_matrix(&M1);

    //let elementwise_product = mz1.mul;
    //let result = elementwise_product;

    //// Run the prover
    non_interactive_sumcheck_prove::<PlutoBaseField>(&poly);

    //// Run the verifier
    //verify(proof);

    println!();
    println!("\x1b[32mSparta(0) Execution Complete\x1b[0m");
}
