use crate::air::{air_index_minus_1, Atmosphere};
use crate::{flat, spherical, Path, PathStepper, RayState, RayStateDerivative};

/// The shape of the simulated Earth
#[derive(Clone, Copy)]
pub enum EarthShape {
    Spherical { radius: f64 },
    Flat,
}

/// Structure storing the shape of the underlying world and the atmospheric model.
#[derive(Clone)]
pub struct Environment {
    pub shape: EarthShape,
    pub atmosphere: Atmosphere,
}

impl Environment {
    /// Returns the refractive index of the air at the given altitude minus 1
    pub fn n_minus_1(&self, h: f64) -> f64 {
        let pressure = self.atmosphere.pressure(h);
        let temperature = self.atmosphere.temperature(h);
        let rh = 0.0;

        air_index_minus_1(530e-9, pressure, temperature, rh)
    }

    /// Returns the refractive index of the air at the given altitude.
    #[inline]
    pub fn n(&self, h: f64) -> f64 {
        self.n_minus_1(h) + 1.0
    }

    /// Returns the derivative of the refractive index of the air with respect to the altitude, at
    /// the given altitude
    #[inline]
    pub fn dn(&self, h: f64) -> f64 {
        let epsilon = 0.01;
        let n1 = self.n_minus_1(h - epsilon);
        let n2 = self.n_minus_1(h + epsilon);
        (n2 - n1) / (2.0 * epsilon)
    }

    /// Returns Some(radius in meters) if the planet model is spherical, or None if it's flat.
    pub fn radius(&self) -> Option<f64> {
        match self.shape {
            EarthShape::Spherical { radius } => Some(radius),
            EarthShape::Flat => None,
        }
    }

    pub(crate) fn calc_derivative_spherical(&self, state: &RayState) -> RayStateDerivative {
        let radius = self.radius().unwrap();
        let dh = state.dh * radius;
        let h = state.h;

        let nr = self.n(h);
        let dnr = self.dn(h);

        let r = h + radius;
        let d2h = dh * dh * dnr / nr + r * r * dnr / nr + 2.0 * dh * dh / r + r;

        RayStateDerivative {
            dx: 1.0,
            dh: state.dh,
            d2h: d2h / radius / radius,
        }
    }

    pub(crate) fn calc_derivative_flat(&self, state: &RayState) -> RayStateDerivative {
        let dh = state.dh;
        let h = state.h;

        let nr = self.n(h);
        let dnr = self.dn(h);

        let d2h = dnr / nr * (1.0 + dh * dh);

        RayStateDerivative { dx: 1.0, dh, d2h }
    }

    /// Returns an object representing a light path.
    ///
    /// The path is defined by 3 parameters:
    /// * `start_h` - the starting altitude of the path in meters
    /// * `start_ang` - the initial angle in radians between the path and the horizontal plane;
    /// -π/2 is down, 0 is horizontal, π/2 is up
    /// * `straight` - `true` if the path should be a straight line, `false` if it should be a ray
    /// affected by the atmosphere
    pub fn cast_ray<'a>(
        &'a self,
        start_h: f64,
        start_ang: f64,
        straight: bool,
    ) -> Box<Path<'a> + 'a> {
        match (straight, self.shape) {
            (true, EarthShape::Flat) => Box::new(flat::Line::from_h_ang(start_h, start_ang)),
            (true, EarthShape::Spherical { .. }) => {
                Box::new(spherical::Line::from_h_ang(self, start_h, start_ang))
            }
            (false, EarthShape::Flat) => Box::new(flat::Ray::from_h_ang(self, start_h, start_ang)),
            (false, EarthShape::Spherical { .. }) => {
                Box::new(spherical::Ray::from_h_ang(self, start_h, start_ang))
            }
        }
    }

    /// Returns an object representing a light path.
    ///
    /// The path is defined by 3 parameters:
    /// * `start_h` - the starting altitude of the path in meters
    /// * `start_ang` - the initial angle in radians between the path and the horizontal plane;
    /// -π/2 is down, 0 is horizontal, π/2 is up
    /// * `straight` - `true` if the path should be a straight line, `false` if it should be a ray
    /// affected by the atmosphere
    pub fn cast_ray_stepper<'a>(
        &'a self,
        start_h: f64,
        start_ang: f64,
        straight: bool,
    ) -> Box<PathStepper<Item = RayState> + 'a> {
        match (straight, self.shape) {
            (true, EarthShape::Flat) => {
                flat::Line::from_h_ang(start_h, start_ang).into_path_stepper()
            }
            (true, EarthShape::Spherical { .. }) => {
                spherical::Line::from_h_ang(self, start_h, start_ang).into_path_stepper()
            }
            (false, EarthShape::Flat) => {
                flat::Ray::from_h_ang(self, start_h, start_ang).into_path_stepper()
            }
            (false, EarthShape::Spherical { .. }) => {
                spherical::Ray::from_h_ang(self, start_h, start_ang).into_path_stepper()
            }
        }
    }

    /// Returns an object representing a light path.
    ///
    /// Instead of using the initial angle, this method chooses a ray that will hit a given target.
    /// The target is defined as distance and altitude.
    ///
    /// * `start_h` - the initial altitude of the path in meters
    /// * `tgt_h` - the altitude of the target point in meters
    /// * `tgt_dist` - the distance of the target point from the initial point, in meters
    /// * `straight` - `true` if the path should be a straight line, `false` if it should be a ray
    /// affected by the atmosphere
    ///
    /// The ray is calculated by performing a binary search on the initial angle.
    pub fn cast_ray_target<'a>(
        &'a self,
        start_h: f64,
        tgt_h: f64,
        tgt_dist: f64,
        straight: bool,
    ) -> Box<Path<'a> + 'a> {
        if straight {
            match self.shape {
                EarthShape::Flat => {
                    Box::new(flat::Line::from_two_points(start_h, 0.0, tgt_h, tgt_dist))
                }
                EarthShape::Spherical { radius } => Box::new(spherical::Line::from_two_points(
                    self,
                    start_h,
                    0.0,
                    tgt_h,
                    tgt_dist / radius,
                )),
            }
        } else {
            let (mut min_ang, mut max_ang) = (-1.5, 1.5);
            let epsilon = 1e-9;

            while max_ang - min_ang > epsilon {
                let cur_ang = 0.5 * (min_ang + max_ang);
                let ray = self.cast_ray(start_h, cur_ang, straight);
                let h = ray.h_at_dist(tgt_dist);
                if h > tgt_h {
                    max_ang = cur_ang;
                } else {
                    min_ang = cur_ang;
                }
            }

            self.cast_ray(start_h, 0.5 * (min_ang + max_ang), straight)
        }
    }
}
