use std::collections::HashMap;
use std::fmt;

#[derive(Clone)]
pub struct FiniteField {
    value: i32,
    prime: i32,
}

impl FiniteField {
    pub fn new(value: i32, prime: i32) -> Self {
        let mut v = value % prime;
        if v < 0 {
            v += prime;
        }
        Self { value: v, prime }
    }

    pub fn get_value(&self) -> i32 {
        self.value
    }

    pub fn add(&self, other: &Self) -> Self {
        Self::new((self.value + other.value) % self.prime, self.prime)
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self::new((self.value - other.value) % self.prime, self.prime)
    }

    pub fn mul(&self, other: &Self) -> Self {
        Self::new((self.value * other.value) % self.prime, self.prime)
    }

    pub fn div(&self, other: &Self) -> Self {
        Self::new((self.value * other.inv().value) % self.prime, self.prime)
    }

    pub fn inv(&self) -> Self {
        Self::new(mod_inv(self.value, self.prime), self.prime)
    }
}

fn mod_inv(a: i32, m: i32) -> i32 {
    let mut mn = (m, a);
    let mut xy = (0, 1);
    while mn.1 != 0 {
        xy = (xy.1, xy.0 - mn.0 / mn.1 * xy.1);
        mn = (mn.1, mn.0 % mn.1);
    }
    while xy.0 < 0 {
        xy.0 += m;
    }
    xy.0
}

impl fmt::Debug for FiniteField {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

pub struct Polynomial {
    coeffs: Vec<FiniteField>,
}

impl Polynomial {
    pub fn new(coeffs: Vec<FiniteField>) -> Self {
        Self { coeffs }
    }

    pub fn add(&self, other: &Self) -> Self {
        let max_len = self.coeffs.len().max(other.coeffs.len());
        let mut result_coeffs = vec![FiniteField::new(0, self.coeffs[0].prime); max_len];
        for i in 0..max_len {
            result_coeffs[i] = if i < self.coeffs.len() { self.coeffs[i].clone() } else { FiniteField::new(0, self.coeffs[0].prime) };
            if i < other.coeffs.len() {
                result_coeffs[i] = result_coeffs[i].add(&other.coeffs[i]);
            }
        }
        Self::new(result_coeffs)
    }

    pub fn sub(&self, other: &Self) -> Self {
        let max_len = self.coeffs.len().max(other.coeffs.len());
        let mut result_coeffs = vec![FiniteField::new(0, self.coeffs[0].prime); max_len];
        for i in 0..max_len {
            result_coeffs[i] = if i < self.coeffs.len() { self.coeffs[i].clone() } else { FiniteField::new(0, self.coeffs[0].prime) };
            if i < other.coeffs.len() {
                result_coeffs[i] = result_coeffs[i].sub(&other.coeffs[i]);
            }
        }
        Self::new(result_coeffs)
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut result_coeffs = vec![FiniteField::new(0, self.coeffs[0].prime); self.coeffs.len() + other.coeffs.len() - 1];
        for i in 0..self.coeffs.len() {
            for j in 0..other.coeffs.len() {
                result_coeffs[i + j] = result_coeffs[i + j].add(&self.coeffs[i].mul(&other.coeffs[j]));
            }
        }
        Self::new(result_coeffs)
    }

    pub fn evaluate(&self, x: &FiniteField) -> FiniteField {
        let mut result = FiniteField::new(0, self.coeffs[0].prime);
        let mut power = FiniteField::new(1, self.coeffs[0].prime);

        for coeff in &self.coeffs {
            result = result.add(&coeff.mul(&power));
            power = power.mul(x);
        }
        result
    }
}

pub struct PolynomialRing {
    variables: HashMap<String, Polynomial>,
}

impl PolynomialRing {
    pub fn new(prime: i32, variable_count: usize, variable_names: &str) -> Self {
        let mut variables = HashMap::new();
        for (index, name) in variable_names.split(", ").enumerate() {
            let mut coeffs = vec![FiniteField::new(0, prime); variable_count];
            if index < variable_count {
                coeffs[index] = FiniteField::new(1, prime);
            }
            variables.insert(name.to_string(), Polynomial::new(coeffs));
        }
        Self { variables }
    }

    pub fn gens(&self) -> Vec<&Polynomial> {
        let mut polys = Vec::new();
        for name in self.variables.keys() {
            polys.push(self.variables.get(name).unwrap());
        }
        polys
    }
}
