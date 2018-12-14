// calculation of the refractive index of air based on
// https://emtoolbox.nist.gov/wavelength/Documentation.asp#ComparisonCiddorandEdlenEquations
// Uses the modified Edlen equation

use super::p_sv;

const A: f64 = 8342.54;
const B: f64 = 2406147.0;
const C: f64 = 15998.0;
const D: f64 = 96095.43;
const E: f64 = 0.601;
const F: f64 = 0.00972;
const G: f64 = 0.003661;

// takes wavelength in meters
#[inline]
fn s(lambda: f64) -> f64 {
    let lambda_um = lambda * 1e6;
    1.0 / lambda_um / lambda_um
}

#[inline]
fn ns_minus_1(lambda: f64) -> f64 {
    let s = s(lambda);
    (A + B / (130.0 - s) + C / (38.9 - s)) * 1e-8
}

#[inline]
fn x(p: f64, t: f64) -> f64 {
    let t = t - 273.15;
    (1.0 + (E - F * t) * p * 1e-8) / (1.0 + G * t)
}

#[inline]
fn ntp_minus_1(p: f64, t: f64, lambda: f64) -> f64 {
    p * ns_minus_1(lambda) * x(p, t) / D
}

/// Returns the air refractive index for the given wavelength (`lambda`), at the given pressure
/// (`p`), temperature (`t`) and relative humidity (`rh`), decreased by 1
pub fn air_index_minus_1(lambda: f64, p: f64, t: f64, rh: f64) -> f64 {
    let pv = rh / 100.0 * p_sv(t);
    ntp_minus_1(p, t, lambda) - (292.75 / t) * (3.7345 - s(lambda) * 0.0401) * pv * 1e-10
}

/// Returns the air refractive index for the given wavelength (`lambda`), at the given pressure
/// (`p`), temperature (`t`) and relative humidity (`rh`)
pub fn air_index(lambda: f64, p: f64, t: f64, rh: f64) -> f64 {
    air_index_minus_1(lambda, p, t, rh) + 1.0
}
