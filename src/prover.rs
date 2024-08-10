extern crate nalgebra as na;
extern crate rand;

use rand::Rng;
use na::{DMatrix, DVector};

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

pub fn row(x1: i32, x2: i32) -> i32 {
    x1 + 2 * x2
}

pub fn col(y1: i32, y2: i32, y3: i32) -> i32 {
    y1 + 2 * y2 + 4 * y3
}

// Equality functions returning 1 if equal, 0 otherwise
pub fn eqx(x1: i32, x2: i32, x11: i32, x22: i32) -> i32 {
    (x1 == x11 && x2 == x22) as i32
}

pub fn eqy(y1: i32, y11: i32, y2: i32, y22: i32, y3: i32, y33: i32) -> i32 {
    (y1 == y11 && y2 == y22 && y3 == y33) as i32
}

pub fn mi_linear(mi: &DMatrix<i32>, x11: i32, x22: i32, y11: i32, y22: i32, y33: i32) -> i32 {
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

pub fn z_linear(zi: &DVector<i32>, y11: i32, y22: i32, y33: i32) -> i32 {
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

pub fn mi_z_prod(mi: &DMatrix<i32>, zi: &DVector<i32>) -> i32 {
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

pub fn compute_g(z1: &DVector<i32>) -> i32 {
    let prod_m1_z1 = mi_z_prod(&M1, z1);
    let prod_m2_z1 = mi_z_prod(&M2, z1);
    let prod_m3_z1 = mi_z_prod(&M3, z1);

    prod_m1_z1 * prod_m2_z1 - prod_m3_z1
}

pub fn compute_h(z1: &DVector<i32>) -> i32 {
    let g = compute_g(z1);
    let mut h = 0;

    for x11 in 0..=1 {
        for x22 in 0..=1 {
            h += g * eqx(x11, x22, x11, x22);
        }
    }

    h
}

pub fn ti_generator(mi: &DMatrix<i32>, zi: &DVector<i32>, r1: i32, r2: i32) -> i32 {
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
pub fn compute_q1(y11: i32, r1: i32, r2: i32, alpha: i32, z1: &DVector<i32>) -> i32 {
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
pub fn compute_q2(y22: i32, r1: i32, r2: i32, r11: i32, alpha: i32, z1: &DVector<i32>) -> i32 {
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
pub fn compute_q3(y33: i32, r1: i32, r2: i32, r11: i32, r22: i32, alpha: i32, z1: &DVector<i32>) -> i32 {
    mi_linear(&M1, r1, r2, r11, r22, y33) * z_linear(z1, r11, r22, y33)
        + alpha * mi_linear(&M2, r1, r2, r11, r22, y33) * z_linear(z1, r11, r22, y33)
        + alpha * alpha * mi_linear(&M3, r1, r2, r11, r22, y33) * z_linear(z1, r11, r22, y33)
}


pub fn prove(z1: DVector<i32>) -> Proof {

    println!();
    println!("\x1b[32mProver Running\x1b[0m");

    // Generate random elements
    let mut rng = rand::thread_rng();

    let r1 = 1;
    let r2 = 0;

    // Create a mutable `Proof` instance
    // TODO refactor this to real code lol
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

    // Transform CCS relation into -> sumcheck
    let prod_m1_z1 = mi_z_prod(&M1, &z1);
    let prod_m2_z1 = mi_z_prod(&M2, &z1);
    let prod_m3_z1 = mi_z_prod(&M3, &z1);

    let g = prod_m1_z1 * prod_m2_z1 - prod_m3_z1;
    let _h = compute_h(&z1);

    let q = g * eqx(proof.beta1, proof.beta2, prod_m1_z1, prod_m2_z1);
    proof.q = q;

    // Generate T values
    let t1 = ti_generator(&M1, &z1, r1, r2);
    let t2 = ti_generator(&M2, &z1, r1, r2);
    let t3 = ti_generator(&M3, &z1, r1, r2);

    proof.t1 = t1;
    proof.t2 = t2;
    proof.t3 = t3;

    let alpha_squared = proof.alpha * proof.alpha;
    proof.t = t1 + proof.alpha * t2 + alpha_squared * t3;

    // Round 1
    let q1: i32 = (0..=1)
        .map(|x2| q * eqx(x2, x2, x2, x2))  // Assuming q depends on x2
        .sum();

    proof.s1 = q1;

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
