# Sparta(0) 

![SPARTA0](./assets/Cover.png)

Sparta(0) is a blazingly fast Rust implementation of the SuperSparta polynomial IOP based on the excellent [SuperSparta by Hand](https://anoma.net/blog/superspartan-by-hand) blog posts.

## Protocol

This code enables the user to convert an arbitrary numerical relation to Customizable Constraint System (CCS), which is a generalized constraint system designed to 
simultaneously capture R1CS, Plonkish, and AIR constraints.

$$
\sum_{i=1}^{q-1}c_i \cdot \bigcirc_{j \in S_i} M_j \cdot z = 0
$$

Where $z$:

$$
z = (w, 1,x) \in \mathbb{F}^n
$$

A relation can be represented in CCS form in the following way: 

$$
c_1 \cdot (M_1 \cdot z \circ M_2 \cdot z) + c_2 \cdot (M_3 \cdot z) = 0
$$

Where $z$ contains the elements of our relation.

We then conduct the following sumcheck reduction on the relation: 

$$
g(a) := \tilde(eq)(\tau, a) \cdot \sum_{i=0}^{q-1} c_i \prod_{j \in S_i} \left( \sum_{y \in {0,1}^{log \ m}} \tilde{M}_j(a, y) \cdot \tilde{Z}(y) \right)
$$ 
    
and confirm that $g$ is equal to 0:
    
$$
\sum_{\mathcal{b} \in \{0,1\}^{log \ m}} g(b) 
$$

## Getting Started

To build and run the repository.

```
cargo run --release
```

## Reference
1. [Superspartan by Hand](https://anoma.net/blog/superspartan-by-hand)
2. [Hypernova by Hand](https://anoma.net/blog/hypernova-by-hand)
3. [Customizable constraint systems for succinct arguments](https://eprint.iacr.org/2023/552.pdf)
4. [HyperNova: Recursive arguments for
customizable constraint systems](https://eprint.iacr.org/2023/573.pdf)
5. [A Time-Space Tradeoff for the Sumcheck Prover](https://eprint.iacr.org/2024/524)






