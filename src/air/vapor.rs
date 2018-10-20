//! Calculation of saturated vapor pressure

const K1: f64 = 1.16705214528e3;
const K2: f64 = -7.24213167032e5;
const K3: f64 = -1.70738469401e1;
const K4: f64 = 1.20208247025e4;
const K5: f64 = -3.23255503223e6;
const K6: f64 = 1.49151086135e1;
const K7: f64 = -4.82326573616e3;
const K8: f64 = 4.05113405421e5;
const K9: f64 = -2.38555575678e-1;
const K10: f64 = 6.50175348448e2;

#[inline]
fn omega(t: f64) -> f64 {
    t + K9 / (t - K10)
}

#[inline]
fn a(t: f64) -> f64 {
    let o = omega(t);
    o * o + K1 * o + K2
}

#[inline]
fn b(t: f64) -> f64 {
    let o = omega(t);
    K3 * o * o + K4 * o + K5
}

#[inline]
fn c(t: f64) -> f64 {
    let o = omega(t);
    K6 * o * o + K7 * o + K8
}

#[inline]
fn x(t: f64) -> f64 {
    let a = a(t);
    let b = b(t);
    let c = c(t);
    -b + (b * b - 4.0 * a * c).sqrt()
}

/// calculates the saturated vapor pressure over water
pub fn p_sv(temp: f64) -> f64 {
    (2.0 * c(temp) / x(temp)).powi(4) * 1e6
}
