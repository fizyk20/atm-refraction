use super::{Path, PathStepper};
use crate::{Environment, RayState};
use na::integration::{Integrator, RK4Integrator, StepSize};

pub struct Line<'a> {
    env: &'a Environment,
    rmin: f64,
    phimin: f64,
}

impl<'a> Line<'a> {
    pub fn from_h_ang(env: &Environment, h: f64, ang: f64) -> Line {
        Line {
            env,
            rmin: (h + env.radius().unwrap()) * ang.cos(),
            phimin: -ang,
        }
    }

    pub fn from_two_points(
        env: &'a Environment,
        h1: f64,
        phi1: f64,
        h2: f64,
        phi2: f64,
    ) -> Line<'a> {
        let r1 = h1 + env.radius().unwrap();
        let r2 = h2 + env.radius().unwrap();
        let a = r1 / r2;
        let tanphi = (a * phi1.cos() - phi2.cos()) / (phi2.sin() - a * phi1.sin());
        let phimin = tanphi.atan();
        Line {
            env,
            rmin: r1 * (phi1 - phimin).cos(),
            phimin,
        }
    }

    pub fn r(&self, phi: f64) -> f64 {
        self.rmin / (phi - self.phimin).cos()
    }
}

impl<'a, 'b: 'a> Path<'a> for Line<'b> {
    fn h_at_dist(&self, dist: f64) -> f64 {
        let r = self.env.radius().unwrap();
        self.r(dist / r) - r
    }

    fn angle_at_dist(&self, dist: f64) -> f64 {
        dist / self.env.radius().unwrap() - self.phimin
    }

    fn into_path_stepper(self) -> Box<dyn PathStepper<Item = RayState> + 'a> {
        Box::new(LineStepper::new(self.env, self, 1.0))
    }
}

pub struct LineStepper<'a> {
    env: &'a Environment,
    x: f64,
    line: Line<'a>,
    step: f64,
}

impl<'a> LineStepper<'a> {
    fn new(env: &'a Environment, line: Line<'a>, step: f64) -> Self {
        Self {
            env,
            x: 0.0,
            line,
            step,
        }
    }

    fn as_state(&self) -> RayState {
        let h = self.line.h_at_dist(self.x);
        let r = self.env.radius().unwrap();
        RayState {
            x: self.x,
            h,
            dh: self.line.angle_at_dist(self.x).tan() * (h + r) / r,
        }
    }
}

impl Iterator for LineStepper<'_> {
    type Item = RayState;

    fn next(&mut self) -> Option<RayState> {
        self.x += self.step;
        Some(self.as_state())
    }
}

impl PathStepper for LineStepper<'_> {
    fn set_step_size(&mut self, step: f64) {
        self.step = step;
    }
}

pub struct Ray<'a> {
    env: &'a Environment,
    start_h: f64,
    start_dh: f64,
}

impl Ray<'_> {
    pub fn from_h_ang(env: &Environment, h: f64, ang: f64) -> Ray {
        let r = env.radius().unwrap();
        let dh = (h + r) * ang.tan() / r;
        Ray {
            env,
            start_h: h,
            start_dh: dh,
        }
    }

    fn state_at_dist(&self, dist: f64) -> RayState {
        let tgt_dist = dist.abs();
        let mut state = RayState {
            x: 0.0,
            h: self.start_h,
            dh: if dist >= 0.0 {
                self.start_dh
            } else {
                -self.start_dh
            },
        };

        let def_step = 5.0;
        let mut integrator = RK4Integrator::new(def_step);
        while state.x < tgt_dist - def_step {
            integrator.propagate_in_place(
                &mut state,
                |state| self.env.calc_derivative_spherical(state),
                StepSize::UseDefault,
            );
        }
        let last_step = tgt_dist - state.x;
        integrator.propagate_in_place(
            &mut state,
            |state| self.env.calc_derivative_spherical(state),
            StepSize::Step(last_step),
        );

        state
    }
}

impl<'a> Path<'a> for Ray<'a> {
    fn h_at_dist(&self, dist: f64) -> f64 {
        let state = self.state_at_dist(dist);
        state.h
    }

    fn angle_at_dist(&self, dist: f64) -> f64 {
        let state = self.state_at_dist(dist);
        state.get_angle(self.env)
    }

    fn into_path_stepper(self) -> Box<dyn PathStepper<Item = RayState> + 'a> {
        let state = RayState {
            x: 0.0,
            h: self.start_h,
            dh: self.start_dh,
        };
        Box::new(RayStepper::new(state, self.env, 1.0))
    }
}

pub struct RayStepper<'a> {
    cur_state: RayState,
    env: &'a Environment,
    integrator: RK4Integrator,
}

impl<'a> RayStepper<'a> {
    fn new(state: RayState, env: &'a Environment, step_size: f64) -> Self {
        Self {
            cur_state: state,
            env,
            integrator: RK4Integrator::new(step_size),
        }
    }
}

impl Iterator for RayStepper<'_> {
    type Item = RayState;

    fn next(&mut self) -> Option<Self::Item> {
        let env = self.env;
        self.integrator.propagate_in_place(
            &mut self.cur_state,
            |state| env.calc_derivative_spherical(state),
            StepSize::UseDefault,
        );
        Some(self.cur_state)
    }
}

impl PathStepper for RayStepper<'_> {
    fn set_step_size(&mut self, step: f64) {
        self.integrator.set_default_step(step);
    }
}
