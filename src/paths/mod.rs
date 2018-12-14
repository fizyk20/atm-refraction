pub mod flat;
pub mod spherical;

use crate::RayState;

pub trait Path<'a> {
    fn h_at_dist(&self, dist: f64) -> f64;
    fn angle_at_dist(&self, dist: f64) -> f64;
    fn into_path_stepper(self) -> Box<PathStepper<Item = RayState> + 'a>;
}

pub trait PathStepper: Iterator {
    fn set_step_size(&mut self, step: f64);
}
