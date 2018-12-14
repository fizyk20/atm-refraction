// calculation of the refractive index of air based on
// https://emtoolbox.nist.gov/wavelength/Documentation.asp#ComparisonCiddorandEdlenEquations

//! A module providing the tooling for atmospheric models.

mod atmosphere;
mod parser;
mod refractive;
mod vapor;

pub use self::atmosphere::{get_atmosphere, us76_atmosphere, Atmosphere};
pub use self::refractive::{air_index, air_index_minus_1};
pub use self::vapor::p_sv;
