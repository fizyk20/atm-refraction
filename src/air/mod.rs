// calculation of the refractive index of air based on
// https://emtoolbox.nist.gov/wavelength/Documentation.asp#ComparisonCiddorandEdlenEquations

mod refractive;
mod vapor;

pub use self::refractive::air_index;
pub use self::vapor::p_sv;
