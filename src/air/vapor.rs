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
fn d_omega(t: f64) -> f64 {
    1.0 - K9 / (t - K10) / (t - K10)
}

#[inline]
fn a(t: f64) -> f64 {
    let o = omega(t);
    o * o + K1 * o + K2
}

#[inline]
fn da(t: f64) -> f64 {
    let o = omega(t);
    let d_o = d_omega(t);
    2.0 * o * d_o + K1 * d_o
}

#[inline]
fn b(t: f64) -> f64 {
    let o = omega(t);
    K3 * o * o + K4 * o + K5
}

#[inline]
fn db(t: f64) -> f64 {
    let o = omega(t);
    let d_o = d_omega(t);
    2.0 * K3 * o * d_o + K4 * d_o
}

#[inline]
fn c(t: f64) -> f64 {
    let o = omega(t);
    K6 * o * o + K7 * o + K8
}

#[inline]
fn dc(t: f64) -> f64 {
    let o = omega(t);
    let d_o = d_omega(t);
    2.0 * K6 * o * d_o + K7 * d_o
}

#[inline]
fn x(t: f64) -> f64 {
    let a = a(t);
    let b = b(t);
    let c = c(t);
    -b + (b * b - 4.0 * a * c).sqrt()
}

#[inline]
fn dx(t: f64) -> f64 {
    let a = a(t);
    let b = b(t);
    let c = c(t);
    let da = da(t);
    let db = db(t);
    let dc = dc(t);
    let delta = b * b - 4.0 * a * c;
    -db + 0.5 / delta.sqrt() * (2.0 * b * db - 4.0 * a * dc - 4.0 * c * da)
}

/// calculates the saturated vapor pressure over water
pub fn p_sv(temp: f64) -> f64 {
    (2.0 * c(temp) / x(temp)).powi(4) * 1e6
}

/// calculates the derivative of the saturated vapor pressure over water with regard to temperature
pub fn dp_sv(temp: f64) -> f64 {
    let c = c(temp);
    let x = x(temp);
    let dc = dc(temp);
    let dx = dx(temp);
    4.0 * (2.0 * c / x).powi(3) * 1e6 * (2.0 * dc / x - 2.0 * c / x / x * dx)
}
