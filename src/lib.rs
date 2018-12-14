//! This crate supplies tooling for calculations of the paths of light rays in the atmosphere.
//!
//! It supports different theoretical shapes of the planet (arbitrary radius, or even flat) and
//! arbitrary atmospheric models (defined by reference temperature and pressure at some altitudes
//! and temperature gradients).
extern crate numeric_algs as na;

/// Module containing tools for defining non-standard atmospheric models.
pub mod air;
mod environment;
mod paths;
mod ray_state;

pub use crate::environment::*;
pub use crate::paths::*;
pub use crate::ray_state::*;
