use crate::air::{air_index_minus_1, Atmosphere};
use crate::{flat, spherical, Path, RayState, RayStateDerivative};

/// the shape of the simulated Earth
#[derive(Clone, Copy)]
pub enum EarthShape {
    Spherical { radius: f64 },
    Flat,
}

#[derive(Clone)]
pub struct Environment {
    pub shape: EarthShape,
    pub atmosphere: Atmosphere,
}

impl Environment {
    pub fn n_minus_1(&self, h: f64) -> f64 {
        let pressure = self.atmosphere.pressure(h);
        let temperature = self.atmosphere.temperature(h);
        let rh = 0.0;

        air_index_minus_1(530e-9, pressure, temperature, rh)
    }

    #[inline]
    pub fn n(&self, h: f64) -> f64 {
        self.n_minus_1(h) + 1.0
    }

    #[inline]
    pub fn dn(&self, h: f64) -> f64 {
        let epsilon = 0.01;
        let n1 = self.n_minus_1(h - epsilon);
        let n2 = self.n_minus_1(h + epsilon);
        (n2 - n1) / (2.0 * epsilon)
    }

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
