use air::Atmosphere;
use na::integration::{Integrator, RK4Integrator, StepSize};
use path::Path;
use ray_state::*;

#[derive(Clone)]
pub struct Ray {
    start_h: f64,
    start_dh: f64,
    atm: Atmosphere,
}

impl Ray {
    pub fn from_h_ang(atm: Atmosphere, h: f64, ang: f64) -> Ray {
        let dh = ang.tan();
        Ray {
            start_h: h,
            start_dh: dh,
            atm,
        }
    }

    fn state_at_dist(&self, dist: f64) -> RayState {
        let tgt_x = if dist >= 0.0 { dist * 1e3 } else { -dist * 1e3 };
        let mut state = RayState {
            x: 0.0,
            h: self.start_h,
            dr: if dist >= 0.0 {
                self.start_dh
            } else {
                -self.start_dh
            },
        };

        let mut integrator = RK4Integrator::new(1.0);
        while state.x < tgt_x {
            integrator.propagate_in_place(
                &mut state,
                |state| calc_derivative_flat(&self.atm, state),
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
        self.start_dh.atan()
    }

    fn h_at_dist(&self, dist: f64) -> f64 {
        let state = self.state_at_dist(dist);
        state.h
    }

    fn angle_at_dist(&self, dist: f64) -> f64 {
        let state = self.state_at_dist(dist);
        state.dr.atan()
    }
}
