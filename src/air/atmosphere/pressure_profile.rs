use std::collections::BTreeMap;

use super::{
    vertical_profile::{VerticalFunction, VerticalProfile},
    A,
};

use cubic_splines::Factors;

#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub enum PressureFunction {
    /// p0 * exp(lambda * (h-h0))
    Exponential { p0: f64, h0: f64, lambda: f64 },
    /// p0 * (1 + a * (h - h0)) ^ exp
    Power { p0: f64, h0: f64, a: f64, exp: f64 },
    /// p0 * (1 + a1 * (h - h0)) ^ exp1 * (1 + a2 * (h - h0)) ^ exp2 * (1 + a3 * (h - h0)) ^ exp3
    /// Used when the temperature is of the form a(h-h1)(h-h2)(h-h3)
    TriplePower {
        p0: f64,
        h0: f64,
        a: [f64; 3],
        exp: [f64; 3],
    },
    /// p0 * (1 + a1 * (h - h0))^exp1 * (a2 * (h - h0)^2 + b2*(h - h0) + 1)^exp2 *
    /// exp(lambda*atan((h - h0)/(a3 * h0 + b3)))
    /// Used when the temperature is of the form a(h-h1)(h^2 + bh + c) (irreducible)
    PowerWithAtan {
        p0: f64,
        h0: f64,
        a1: f64,
        exp1: f64,
        a2: f64,
        b2: f64,
        exp2: f64,
        lambda: f64,
        a3: f64,
        b3: f64,
    },
}

impl PressureFunction {
    pub fn eval(&self, h: f64) -> f64 {
        match *self {
            PressureFunction::Exponential { p0, h0, lambda } => p0 * (lambda * (h - h0)).exp(),
            PressureFunction::Power { p0, h0, a, exp } => p0 * (1.0 + a * (h - h0)).powf(exp),
            PressureFunction::TriplePower { p0, h0, a, exp } => {
                p0 * (1.0 + a[0] * (h - h0)).powf(exp[0])
                    * (1.0 + a[1] * (h - h0)).powf(exp[1])
                    * (1.0 + a[2] * (h - h0)).powf(exp[2])
            }
            PressureFunction::PowerWithAtan {
                p0,
                h0,
                a1,
                exp1,
                a2,
                b2,
                exp2,
                lambda,
                a3,
                b3,
            } => {
                p0 * (1.0 + a1 * (h - h0)).powf(exp1)
                    * (1.0 + b2 * (h - h0) + a2 * (h - h0) * (h - h0)).powf(exp2)
                    * (lambda * ((h - h0) / (a3 * (h - h0) + b3)).atan()).exp()
            }
        }
    }

    pub fn from_temperature_function(temp_function: &VerticalFunction, p0: f64, h0: f64) -> Self {
        match *temp_function {
            VerticalFunction::Linear { a, b } => {
                if a == 0.0 {
                    PressureFunction::Exponential {
                        p0,
                        h0,
                        lambda: -A / b,
                    }
                } else {
                    PressureFunction::Power {
                        p0,
                        h0,
                        a: a / (a * h0 + b),
                        exp: -A / a,
                    }
                }
            }
            VerticalFunction::Cubic(poly) => match poly.factors() {
                Factors::ThreeLinear {
                    a,
                    x1: h1,
                    x2: h2,
                    x3: h3,
                } => {
                    let v = [
                        1.0 / (h1 - h2) / (h1 - h3),
                        1.0 / (h2 - h1) / (h2 - h3),
                        1.0 / (h3 - h1) / (h3 - h2),
                    ];
                    let exp = [-A * v[0] / a, -A * v[1] / a, -A * v[2] / a];
                    let a = [1.0 / (h0 - h1), 1.0 / (h0 - h2), 1.0 / (h0 - h3)];
                    PressureFunction::TriplePower { p0, h0, a, exp }
                }
                Factors::LinearAndQuadratic { a, x1: h1, b, c } => {
                    let u = h1 * h1 + b * h1 + c;
                    let v = [1.0 / u, -1.0 / u, -(h1 + b) / u];
                    let a1 = 1.0 / (h0 - h1);
                    let exp1 = -A * v[0] / a;
                    let a2 = 1.0 / (h0 * h0 + b * h0 + c);
                    let two_h_b = 2.0 * h0 + b;
                    let b2 = two_h_b * a2;
                    let exp2 = -A * v[1] / 2.0 / a;
                    let sqrt = (4.0 * c - b * b).sqrt();
                    let lambda = -A * (2.0 * v[2] - v[1] * b) / a / sqrt;
                    let a3 = two_h_b / sqrt;
                    let b3 = two_h_b * two_h_b / 2.0 / sqrt + sqrt / 2.0;
                    PressureFunction::PowerWithAtan {
                        p0,
                        h0,
                        a1,
                        exp1,
                        a2,
                        b2,
                        exp2,
                        lambda,
                        a3,
                        b3,
                    }
                }
            },
        }
    }
}

#[derive(Clone, Debug)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct PressureProfile {
    altitude_interval_ends: Vec<f64>,
    pressure_functions: Vec<PressureFunction>,
}

impl PressureProfile {
    pub fn from_temperature_profile(temp: &VerticalProfile, p0: f64, h0: f64) -> Self {
        let (altitude_interval_ends, interval_functions) = temp.internals();
        let (start_index, mut map) = match altitude_interval_ends
            .binary_search_by(|h| h.partial_cmp(&h0).unwrap())
        {
            Ok(index) | Err(index) => {
                let function =
                    PressureFunction::from_temperature_function(&interval_functions[index], p0, h0);
                let mut map = BTreeMap::new();
                let _ = map.insert(index, function);
                (index, map)
            }
        };
        if let Some(start_index_below) = start_index.checked_sub(1) {
            for index in (0..=start_index_below).rev() {
                let h0 = altitude_interval_ends[index];
                let p0 = map[&(index + 1)].eval(h0);
                let _ = map.insert(
                    index,
                    PressureFunction::from_temperature_function(&interval_functions[index], p0, h0),
                );
            }
        }
        if let Some(start_index_above) =
            (start_index < altitude_interval_ends.len()).then(|| start_index + 1)
        {
            for index in start_index_above..interval_functions.len() {
                let h0 = altitude_interval_ends[index - 1];
                let p0 = map[&(index - 1)].eval(h0);
                let _ = map.insert(
                    index,
                    PressureFunction::from_temperature_function(&interval_functions[index], p0, h0),
                );
            }
        }

        let pressure_functions = map.into_iter().map(|(_, fun)| fun).collect();

        PressureProfile {
            altitude_interval_ends: altitude_interval_ends.clone(),
            pressure_functions,
        }
    }

    pub fn eval(&self, h: f64) -> f64 {
        match self
            .altitude_interval_ends
            .binary_search_by(|a| a.partial_cmp(&h).unwrap())
        {
            Ok(index) | Err(index) => self.pressure_functions[index].eval(h),
        }
    }
}
