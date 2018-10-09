use na::integration::{Integrator, RK4Integrator, StepSize};
use path::Path;
use ray_state::*;

pub struct Ray {
    radius: f64,
    start_h: f64,
    start_dh: f64,
}

impl Ray {
    pub fn from_h_ang(radius: f64, h: f64, ang: f64) -> Ray {
        let dh = (h + radius) * ang.tan();
        Ray {
            radius,
            start_h: h,
            start_dh: dh,
        }
    }

    fn state_at_dist(&self, dist: f64) -> RayState {
        let tgt_phi = if dist >= 0.0 {
            dist * 1e3 / self.radius
        } else {
            -dist * 1e3 / self.radius
        };
        let mut state = RayState {
            x: 0.0,
            h: self.start_h,
            dr: if dist >= 0.0 {
                self.start_dh
            } else {
                -self.start_dh
            },
        };

        let mut integrator = RK4Integrator::new(1e-5);
        while state.x < tgt_phi {
            integrator.propagate_in_place(
                &mut state,
                |state| calc_derivative_spherical(self.radius, state),
                StepSize::UseDefault,
            );
        }

        state
    }
}

impl Path for Ray {
    fn start_h(&self) -> f64 {
        self.start_h
    }

    fn start_angle(&self) -> f64 {
        (self.start_dh / (self.start_h + self.radius)).atan()
    }

    fn h_at_dist(&self, dist: f64) -> f64 {
        let state = self.state_at_dist(dist);
        state.h
    }

    fn angle_at_dist(&self, dist: f64) -> f64 {
        let state = self.state_at_dist(dist);
        (state.dr / (state.h + self.radius)).atan()
    }
}
