#[macro_use]
extern crate lazy_static;
extern crate nalgebra as na;
extern crate rand;
extern crate num_bigint;
extern crate num_traits;

mod polynomial;

use polynomial::{FiniteField, PolynomialRing};
use na::{DMatrix, DVector};
use rand::Rng;

pub struct Proof {
    pub alpha: i32,
    pub beta1: i32,
    pub beta2: i32,
    pub r1: i32,
    pub r2: i32,
    pub r11: i32,
    pub r22: i32,
    pub r33: i32,
    pub z: DVector<i32>,
    pub q: i32,
    pub s1: i32,
    pub s2: i32,
    pub s3: i32,
    pub t: i32,
    pub t1: i32,
    pub t2: i32,
    pub t3: i32,
    pub c: i32,
}

lazy_static! {
    pub static ref M1: DMatrix<i32> = DMatrix::from_row_slice(4, 8, &[
        1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
    ]);

    pub static ref M2: DMatrix<i32> = DMatrix::from_row_slice(4, 8, &[
        0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
    ]);

    pub static ref M3: DMatrix<i32> = DMatrix::from_row_slice(4, 8, &[
        0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
    ]);
}

fn row(x1: i32, x2: i32) -> i32 {
    x1 + 2 * x2
}

fn col(y1: i32, y2: i32, y3: i32) -> i32 {
    y1 + 2 * y2 + 4 * y3
}

// Equality functions returning 1 if equal, 0 otherwise
fn eqx(x1: i32, x2: i32, x11: i32, x22: i32) -> i32 {
    (x1 == x11 && x2 == x22) as i32
}

fn eqy(y1: i32, y11: i32, y2: i32, y22: i32, y3: i32, y33: i32) -> i32 {
    (y1 == y11 && y2 == y22 && y3 == y33) as i32
}

fn mi_linear(mi: &DMatrix<i32>, x11: i32, x22: i32, y11: i32, y22: i32, y33: i32) -> i32 {
    let mut result = 0;

    for x1 in 0..=1 {
        for x2 in 0..=1 {
            for y1 in 0..=1 {
                for y2 in 0..=1 {
                    for y3 in 0..=1 {
                        let row_index = row(x1, x2) as usize;
                        let col_index = col(y1, y2, y3) as usize;
                        result += mi[(row_index, col_index)]
                            * eqx(x1, x2, x11, x22)
                            * eqy(y1, y11, y2, y22, y3, y33);
                    }
                }
            }
        }
    }

    result
}

fn z_linear(zi: &DVector<i32>, y11: i32, y22: i32, y33: i32) -> i32 {
    let mut result = 0;

    for y1 in 0..=1 {
        for y2 in 0..=1 {
            for y3 in 0..=1 {
                let col_index = col(y1, y2, y3) as usize;
                result += zi[col_index] * eqy(y1, y11, y2, y22, y3, y33);
            }
        }
    }

    result
}

fn mi_z_prod(mi: &DMatrix<i32>, zi: &DVector<i32>) -> i32 {
    let mut result = 0;

    for y1 in 0..=1 {
        for y2 in 0..=1 {
            for y3 in 0..=1 {
                result += mi_linear(mi, y1, y2, y1, y2, y3)
                        * z_linear(zi, y1, y2, y3);
            }
        }
    }

    result
}

fn compute_g(z1: &DVector<i32>) -> i32 {
    let prod_m1_z1 = mi_z_prod(&M1, z1);
    let prod_m2_z1 = mi_z_prod(&M2, z1);
    let prod_m3_z1 = mi_z_prod(&M3, z1);

    prod_m1_z1 * prod_m2_z1 - prod_m3_z1
}

fn compute_h(z1: &DVector<i32>) -> i32 {
    let g = compute_g(z1);
    let mut h = 0;

    for x11 in 0..=1 {
        for x22 in 0..=1 {
            h += g * eqx(x11, x22, x11, x22);
        }
    }

    h
}

fn ti_generator(mi: &DMatrix<i32>, zi: &DVector<i32>, r1: i32, r2: i32) -> i32 {
    let mut result = 0;

    for y1 in 0..=1 {
        for y2 in 0..=1 {
            for y3 in 0..=1 {
                let mi_val = mi_linear(mi, r1, r2, y1, y2, y3);
                let zi_val = z_linear(zi, y1, y2, y3);
                result += mi_val * zi_val;
            }
        }
    }

    result
}

// Define q1 as a function that takes y11 as parameter
fn compute_q1(y11: i32, r1: i32, r2: i32, alpha: i32, z1: &DVector<i32>) -> i32 {
    (0..=1)
        .flat_map(|y2| {
            (0..=1).map(move |y3| {
                let mi1 = mi_linear(&M1, r1, r2, y11, y2, y3);
                let mi2 = mi_linear(&M2, r1, r2, y11, y2, y3);
                let mi3 = mi_linear(&M3, r1, r2, y11, y2, y3);
                let z1_val = z_linear(&z1, y11, y2, y3);

                mi1 * z1_val + alpha * mi2 * z1_val + alpha * alpha * mi3 * z1_val
            })
        })
        .sum()
}

// Define compute_q2 as a function that takes y22 as a parameter
fn compute_q2(y22: i32, r1: i32, r2: i32, r11: i32, alpha: i32, z1: &DVector<i32>) -> i32 {
    (0..=1)
        .map(|y3| {
            let mi1 = mi_linear(&M1, r1, r2, r11, y22, y3);
            let mi2 = mi_linear(&M2, r1, r2, r11, y22, y3);
            let mi3 = mi_linear(&M3, r1, r2, r11, y22, y3);
            let z1_val = z_linear(z1, r11, y22, y3);

            mi1 * z1_val + alpha * mi2 * z1_val + alpha * alpha * mi3 * z1_val
        })
        .sum()
}

// Define compute_q3 as a function that takes y33 as a parameter
fn compute_q3(y33: i32, r1: i32, r2: i32, r11: i32, r22: i32, alpha: i32, z1: &DVector<i32>) -> i32 {
    mi_linear(&M1, r1, r2, r11, r22, y33) * z_linear(z1, r11, r22, y33)
        + alpha * mi_linear(&M2, r1, r2, r11, r22, y33) * z_linear(z1, r11, r22, y33)
        + alpha * alpha * mi_linear(&M3, r1, r2, r11, r22, y33) * z_linear(z1, r11, r22, y33)
}

fn main() {
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

    // Run the prover
    let proof = prove(z1);

    // Run the verifier
    verify(proof);
}

fn prove(z1: DVector<i32>) -> Proof {

    println!();
    println!("\x1b[32mProver Running\x1b[0m");

    // Generate random elements
    let mut rng = rand::thread_rng();

    let r1 = 1;
    let r2 = 0;

    // Create a mutable `Proof` instance
    let mut proof = Proof {
        alpha: rng.gen(),
        beta1: rng.gen(),
        beta2: rng.gen(),
        r1: 1,
        r2: 0,
        r11: rng.gen_range(0..=1),
        r22: rng.gen_range(0..=1),
        r33: rng.gen_range(0..=1),
        z: z1.clone(),
        q: 0,
        s1: 0,
        s2: 0,
        s3: 0,
        t: 0,
        t1: 0,
        t2: 0,
        t3: 0,
        c: 0,
    };

    let mz1 = M1.clone() * (&z1);
    let mz2 = M2.clone() * (&z1);
    let mz3 = M3.clone() * (&z1);

    let elementwise_product = mz1.component_mul(&mz2);
    let result = elementwise_product - mz3;

    // CCS relation check
    if result == DVector::from_element(4, 0) {
        let m = M1.nrows();
        let n = M1.ncols();
        println!("      CCS relation check passed.");
        println!("      (m, n) = ({}, {})", m, n);
    } else {
        println!("      CCS relation check failed.");
    }

    // Transform CCS relation into -> sum-check
    let mi_result = mi_linear(&M1, 0, 0, 1, 0, 0);
    let z_result = z_linear(&z1, 1, 0, 0);
    let mi_z_result = mi_z_prod(&M1, &z1);

    let prod_m1_z1 = mi_z_prod(&M1, &z1);
    let prod_m2_z1 = mi_z_prod(&M2, &z1);
    let prod_m3_z1 = mi_z_prod(&M3, &z1);

    let g = prod_m1_z1 * prod_m2_z1 - prod_m3_z1;
    let h = compute_h(&z1);

    let q = g * eqx(proof.beta1, proof.beta2, prod_m1_z1, prod_m2_z1);
    proof.q = q;

    // Round 1
    let q1: i32 = (0..=1)
        .map(|x2| q * eqx(x2, x2, x2, x2))  // Assuming q depends on x2
        .sum();

    proof.s1 = q1;

    let t1 = ti_generator(&M1, &z1, r1, r2);
    let t2 = ti_generator(&M2, &z1, r1, r2);
    let t3 = ti_generator(&M3, &z1, r1, r2);

    proof.t1 = t1;
    proof.t2 = t2;
    proof.t3 = t3;

    let alpha_squared = proof.alpha * proof.alpha;
    proof.t = t1 + proof.alpha * t2 + alpha_squared * t3;

    // Round 2
    let q2_y22_0 = compute_q2(0, r1, r2, proof.r11, proof.alpha, &z1);
    let q2_y22_1 = compute_q2(1, r1, r2, proof.r11, proof.alpha, &z1);

    proof.s2 = q2_y22_0 + q2_y22_1;

    // Round 3
    let q3_y33_0 = compute_q3(0, proof.r1, proof.r2, proof.r11, proof.r22, proof.alpha, &z1);
    let q3_y33_1 = compute_q3(1, proof.r1, proof.r2, proof.r11, proof.r22, proof.alpha, &z1);

    proof.s3 = q3_y33_0 + q3_y33_1;

    // Compute c1
    proof.c = mi_linear(&M1, proof.r1, proof.r2, proof.r11, proof.r22, proof.r33) * z_linear(&z1, proof.r11, proof.r22, proof.r33)
        + proof.alpha * mi_linear(&M2, proof.r1, proof.r2, proof.r11, proof.r22, proof.r33) * z_linear(&z1, proof.r11, proof.r22, proof.r33)
        + proof.alpha * proof.alpha * mi_linear(&M3, proof.r1, proof.r2, proof.r11, proof.r22, proof.r33) * z_linear(&z1, proof.r11, proof.r22, proof.r33);

    return proof;
}

fn verify(proof: Proof) {

    println!();
    println!("\x1b[32mVerifier Running\x1b[0m");
    // Outer sum-check
    let sum_q: i32 = (0..=1)
        .map(|x1| {
            (0..=1)
                .map(|x2| proof.q * eqx(x1, x2, x1, x2))  // Assuming q depends on x1 and x2
                .sum::<i32>()
        })
        .sum();

    let result = sum_q == 0;
    println!("      Initial Check: Outer Sum-check: Assertion {:?}", result);

    // Round 1
    let q1_y11_0 = compute_q1(0, proof.r1, proof.r2, proof.alpha, &proof.z);
    let q1_y11_1 = compute_q1(1, proof.r1, proof.r2, proof.alpha, &proof.z);

    let result = q1_y11_0 + q1_y11_1 == proof.t;
    println!("      Round 1: Assertion (q1(y11=0) + q1(y11=1) == T) holds: {}", result);

    // Round 2
    let q1_r11 = compute_q1(proof.r11, proof.r1, proof.r2, proof.alpha, &proof.z);

    let assertion_result = proof.s2 == q1_r11;
    println!("      Round 2: Assertion (q2(y22=0) + q2(y22=1) == q1(y11=r11)) holds: {}", assertion_result);

    // Round 3
    let q2_r22 = compute_q2(proof.r22, proof.r1, proof.r2, proof.r11, proof.alpha, &proof.z);

    let assertion_result = proof.s3 == q2_r22;
    println!("      Round 3: Assertion (q3(y33=0) + q3(y33=1) == q2(y22=r22)) holds: {}", assertion_result);

    // Final Check
    let q3_r33 = compute_q3(proof.r33, proof.r1, proof.r2, proof.r11, proof.r22, proof.alpha, &proof.z);

    // Assert equality
    let assertion_result = proof.c == q3_r33;
    println!("      Final Check: C1 Assertion (c1 == q3(y33=r33)) holds: {}", assertion_result);

    // Compute Q and (T1 * T2 - T3) * eqx
    let eqx_val = eqx(proof.beta1, proof.beta2, proof.r1, proof.r2);
    let t_product = (proof.t1 * proof.t2 - proof.t3) * eqx_val;

    // Final assert equality
    let assertion_result = proof.q == t_product;
    println!("      Final Check: Q Assertion (Q(x11=r1, x22=r2) == (T1 * T2 - T3) * eqx(x11=r1, x22=r2, x1=beta1, x2=beta2)) holds: {}", assertion_result);
}
