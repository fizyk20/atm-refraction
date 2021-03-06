// calculation of the refractive index of air based on
// https://emtoolbox.nist.gov/wavelength/Documentation.asp#ComparisonCiddorandEdlenEquations

//! A module providing the tooling for atmospheric models.

pub mod atmosphere;
mod refractive;
mod vapor;

pub use self::atmosphere::{us76_atmosphere, Atmosphere, AtmosphereDef};
pub use self::refractive::{air_index, d_air_index};
pub use self::vapor::{dp_sv, p_sv};
