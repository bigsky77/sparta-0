use crate::prover::{Proof, compute_q1, compute_q2, compute_q3, eqx};

pub fn verify(proof: Proof) {

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
