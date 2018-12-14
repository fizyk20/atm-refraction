pub(crate) mod flat;
pub(crate) mod spherical;

use crate::RayState;

/// The trait representing a light path.
pub trait Path<'a> {
    /// Returns the altitude (in meters) at which the path is passing at the given distance (in
    /// meters) from the initial point.
    fn h_at_dist(&self, dist: f64) -> f64;
    /// Returns the angle (in radians) between the path and the horizontal plane at the given
    /// distance (in meters) from the initial point.
    fn angle_at_dist(&self, dist: f64) -> f64;
    /// Returns a "stepper" - an iterator that performs one integration step along the path on
    /// every call to `next()`
    fn into_path_stepper(self) -> Box<PathStepper<Item = RayState> + 'a>;
}

/// The trait representing a "stepper" - an iterator performing one integration step along the
/// path on every call to `next()`
pub trait PathStepper: Iterator {
    /// Sets the step size for the iterations
    fn set_step_size(&mut self, step: f64);
}
