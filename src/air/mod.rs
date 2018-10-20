// calculation of the refractive index of air based on
// https://emtoolbox.nist.gov/wavelength/Documentation.asp#ComparisonCiddorandEdlenEquations

mod refractive;
mod us76;
mod vapor;

pub use self::refractive::{air_index, air_index_minus_1};
pub use self::us76::{pressure, temperature};
pub use self::vapor::p_sv;
