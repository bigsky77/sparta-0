#[macro_use]
extern crate lazy_static;
extern crate nalgebra as na;
extern crate rand;

use na::{DMatrix, DVector};
use rand::Rng;

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
    lazy_static::lazy_static! {
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
    // Generate random elements beta1 and beta2
    let mut rng = rand::thread_rng();
    let beta1: i32 = rng.gen();
    let beta2: i32 = rng.gen();

    println!("Random Elements: (beta1, beta2) = ({}, {})", beta1, beta2);

    let x = vec![0, 1, 1, 2, 3, 6, 6];
    let w: Vec<i32> = Vec::new();
    let mut combined = x.clone();
    combined.extend(&w);
    combined.push(1);
    let z1 = DVector::from_vec(combined);

    let mz1 = M1.clone() * (&z1);
    let mz2 = M2.clone() * (&z1);
    let mz3 = M3.clone() * (&z1);

    let elementwise_product = mz1.component_mul(&mz2);
    let result = elementwise_product - mz3;

    if result == DVector::from_element(4, 0) {
        let m = M1.nrows();
        let n = M1.ncols();
        println!("CCS relation check passed.");
        println!("(m, n) = ({}, {})", m, n);
    } else {
        println!("CCS relation check failed.");
    }

    let mi_result = mi_linear(&M1, 0, 0, 1, 0, 0);
    println!("Result: {:?}", mi_result);

    let z_result = z_linear(&z1, 1, 0, 0);
    println!("Result: {:?}", z_result);

    let mi_z_result = mi_z_prod(&M1, &z1);
    println!("Result: {:?}", mi_z_result);

    let prod_m1_z1 = mi_z_prod(&M1, &z1);
    let prod_m2_z1 = mi_z_prod(&M2, &z1);
    let prod_m3_z1 = mi_z_prod(&M3, &z1);

    let g = prod_m1_z1 * prod_m2_z1 - prod_m3_z1;
    println!("G: {:?}", g);

    let h = compute_h(&z1);
    println!("h: {:?}", h);

    let q = g * eqx(beta1, beta2, prod_m1_z1, prod_m2_z1);
    println!("Q: {:?}", q);

    let sum_q: i32 = (0..=1)
        .map(|x1| {
            (0..=1)
                .map(|x2| q * eqx(x1, x2, x1, x2))  // Assuming q depends on x1 and x2
                .sum::<i32>()
        })
        .sum();

    let result = sum_q == 0;
    println!("Result: {:?}", result);

    let q1: i32 = (0..=1)
        .map(|x2| q * eqx(x2, x2, x2, x2))  // Assuming q depends on x2
        .sum();

    let s1 = q1;
    println!("S1: {:?}", s1);

    let r1 = 1;
    let r2 = 0;

    let t1 = ti_generator(&M1, &z1, r1, r2);
    let t2 = ti_generator(&M2, &z1, r1, r2);
    let t3 = ti_generator(&M3, &z1, r1, r2);

    println!("Ti: {:?}", t1);
    println!("Ti: {:?}", t2);
    println!("Ti: {:?}", t3);

    let alpha: i32 = rng.gen();

    let alpha_squared = alpha * alpha;
    let t = t1 + alpha * t2 + alpha_squared * t3;

    println!("T: {:?}", t);

    let mut f1 = 0;

    for y2 in 0..=1 {
        for y3 in 0..=1 {
            let mi1 = mi_linear(&M1, r1, r2, y2, y2, y3);
            let mi2 = mi_linear(&M2, r1, r2, y2, y2, y3);
            let mi3 = mi_linear(&M3, r1, r2, y2, y2, y3);
            let z1_val = z_linear(&z1, y2, y2, y3);

            f1 += mi1 * z1_val + alpha * mi2 * z1_val + alpha * alpha * mi3 * z1_val;
        }
    }

    let q1 = f1;
    println!("Q1: {:?}", q1);

    let q1_y11_0 = compute_q1(0, r1, r2, alpha, &z1);
    let q1_y11_1 = compute_q1(1, r1, r2, alpha, &z1);

    let result = q1_y11_0 + q1_y11_1 == t;
    println!("Condition (q1(y11=0) + q1(y11=1) == T) holds: {}", result);

    let r11: i32 = rng.gen_range(0..=1); // Generates a random value for r11 within the range [0, 1]
    let r22: i32 = rng.gen_range(0..=1); // Generates a random value for r22 within the range [0, 1]
    let r33: i32 = rng.gen_range(0..=1); // Generates a random value for r22 within the range [0, 1]

    // Compute f2
    let f2: i32 = (0..=1)
        .map(|y3| {
            let mi1 = mi_linear(&M1, r1, r2, r11, r11, y3);
            let mi2 = mi_linear(&M2, r1, r2, r11, r11, y3);
            let mi3 = mi_linear(&M3, r1, r2, r11, r11, y3);
            let z1_val = z_linear(&z1, r11, r11, y3);

            mi1 * z1_val + alpha * mi2 * z1_val + alpha * alpha * mi3 * z1_val
        })
        .sum();

    let q2 = f2;
    println!("q2: {}", q2); // Print the result

    let q2_y22_0 = compute_q2(0, r1, r2, r11, alpha, &z1);
    let q2_y22_1 = compute_q2(1, r1, r2, r11, alpha, &z1);

    let q2_sum = q2_y22_0 + q2_y22_1;

    let q1_r11 = compute_q1(r11, r1, r2, alpha, &z1);

    let assertion_result = q2_sum == q1_r11;
    println!("Assertion (q2(y22=0) + q2(y22=1) == q1(y11=r11)) holds: {}", assertion_result);

    println!("q2(y22=0): {}, q2(y22=1): {}", q2_y22_0, q2_y22_1);
    println!("q2_sum: {}", q2_sum);
    println!("q1(y11=r11): {}", q1_r11);

    // Compute f3
    let f3 = mi_linear(&M1, r1, r2, r11, r22, 0) * z_linear(&z1, r11, r22, 0)
        + alpha * mi_linear(&M2, r1, r2, r11, r22, 0) * z_linear(&z1, r11, r22, 0)
        + alpha * alpha * mi_linear(&M3, r1, r2, r11, r22, 0) * z_linear(&z1, r11, r22, 0);

    let q3 = f3;

    let q3_y33_0 = compute_q3(0, r1, r2, r11, r22, alpha, &z1);
    let q3_y33_1 = compute_q3(1, r1, r2, r11, r22, alpha, &z1);

    let q3_sum = q3_y33_0 + q3_y33_1;

    let q2_r22 = compute_q2(r22, r1, r2, r11, alpha, &z1);

    let assertion_result = q3_sum == q2_r22;
    println!("Assertion (q3(y33=0) + q3(y33=1) == q2(y22=r22)) holds: {}", assertion_result);

    // Compute c1
    let c1 = mi_linear(&M1, r1, r2, r11, r22, r33) * z_linear(&z1, r11, r22, r33)
        + alpha * mi_linear(&M2, r1, r2, r11, r22, r33) * z_linear(&z1, r11, r22, r33)
        + alpha * alpha * mi_linear(&M3, r1, r2, r11, r22, r33) * z_linear(&z1, r11, r22, r33);

    // Compute q3(y33=r33)
    let q3_r33 = compute_q3(r33, r1, r2, r11, r22, alpha, &z1);

    // Assert equality
    let assertion_result = c1 == q3_r33;
    println!("Assertion (c1 == q3(y33=r33)) holds: {}", assertion_result);

    // Compute Q and (T1 * T2 - T3) * eqx
    let eqx_val = eqx(beta1, beta2, r1, r2);
    let t_product = (t1 * t2 - t3) * eqx_val;

    // Final assert equality
    let assertion_result = q == t_product;
    println!("Assertion (Q(x11=r1, x22=r2) == (T1 * T2 - T3) * eqx(x11=r1, x22=r2, x1=beta1, x2=beta2)) holds: {}", assertion_result);

    println!("(T1 * T2 - T3) * eqx(x11=r1, x22=r2, x1=beta1, x2=beta2): {}", t_product);
}
