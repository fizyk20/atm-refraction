// calculation of the refractive index of air based on
// https://emtoolbox.nist.gov/wavelength/Documentation.asp#ComparisonCiddorandEdlenEquations
// Uses the modified Edlen equation

use super::{dp_sv, p_sv};

const A: f64 = 8342.54;
const B: f64 = 2406147.0;
const C: f64 = 15998.0;
const D: f64 = 96095.43;
const E: f64 = 0.601;
const F: f64 = 0.00972;
const G: f64 = 0.003661;

/// Returns the air refractive index for the given wavelength (`lambda`), at the given pressure
/// (`p`), temperature (`t`) and relative humidity (`rh`)
pub fn air_index(lambda: f64, p: f64, t: f64, rh: f64) -> f64 {
    let lambda_um = lambda * 1e6;
    let s = 1.0 / lambda_um / lambda_um;
    let t1 = t - 273.15;

    let alpha = 1e-8 * (A + B / (130.0 - s) + C / (38.9 - s));
    let beta = 1e-8 * E;
    let gamma = -1e-8 * F;
    let delta = D;
    let epsilon = D * G;
    let zeta = (3.7345 - s * 0.0401) * 1e-10;

    let pv = rh / 100.0 * p_sv(t);

    1.0 + alpha * p * (1.0 + beta * p + gamma * t1 * p) / (delta + epsilon * t1)
        - (292.75 / t) * zeta * pv
}

/// Returns the derivative of the air refractive index for the given wavelength (`lambda`) as a
/// function of pressure (`p`), temperature (`t`), relative humidity (`rh`) and their derivatives
/// (`dp`, `dt`, `drh`)
pub fn d_air_index(lambda: f64, p: f64, t: f64, rh: f64, dp: f64, dt: f64, drh: f64) -> f64 {
    let lambda_um = lambda * 1e6;
    let s = 1.0 / lambda_um / lambda_um;
    let t1 = t - 273.15;

    let alpha = 1e-8 * (A + B / (130.0 - s) + C / (38.9 - s));
    let beta = 1e-8 * E;
    let gamma = -1e-8 * F;
    let delta = D;
    let epsilon = D * G;
    let zeta = (3.7345 - s * 0.0401) * 1e-10;

    let pv = rh / 100.0 * p_sv(t);
    let dpv = drh / 100.0 * p_sv(t) + rh / 100.0 * dp_sv(t) * dt;

    alpha * dp * (1.0 + beta * p + gamma * t1 * p) / (delta + epsilon * t1)
        + alpha
            * p
            * ((beta * dp + gamma * t1 * dp + gamma * p * dt) * (delta + epsilon * t1)
                - epsilon * dt * (1.0 + beta * p + gamma * t1 * p))
            / (delta + epsilon * t1)
            / (delta + epsilon * t1)
        + 292.75 / t / t * dt * zeta * pv
        - 292.75 / t * zeta * dpv
}
