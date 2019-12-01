use super::{Path, PathStepper};
use crate::{Environment, RayState};
use na::integration::{Integrator, RK4Integrator, StepSize};

pub struct Line {
    a: f64,
    b: f64,
}

impl Line {
    pub fn from_h_ang(h: f64, ang: f64) -> Line {
        let a = ang.tan();
        Line { a, b: h }
    }

    pub fn from_two_points(h1: f64, x1: f64, h2: f64, x2: f64) -> Line {
        let a = (h2 - h1) / (x2 - x1);
        let b = h1 - a * x1;
        Line { a, b }
    }
}

impl<'a> Path<'a> for Line {
    fn h_at_dist(&self, dist: f64) -> f64 {
        self.a * dist + self.b
    }

    fn angle_at_dist(&self, _dist: f64) -> f64 {
        self.a.atan()
    }

    fn into_path_stepper(self) -> Box<dyn PathStepper<Item = RayState> + 'a> {
        Box::new(LineStepper::new(self, 1.0))
    }
}

pub struct LineStepper {
    x: f64,
    line: Line,
    step: f64,
}

impl LineStepper {
    fn new(line: Line, step: f64) -> Self {
        Self { x: 0.0, line, step }
    }

    fn into_state(&self) -> RayState {
        let h = self.line.h_at_dist(self.x);
        RayState {
            x: self.x,
            h,
            dh: self.line.angle_at_dist(self.x).tan(),
        }
    }
}

impl Iterator for LineStepper {
    type Item = RayState;

    fn next(&mut self) -> Option<RayState> {
        self.x += self.step;
        Some(self.into_state())
    }
}

impl PathStepper for LineStepper {
    fn set_step_size(&mut self, step: f64) {
        self.step = step;
    }
}

pub struct Ray<'a> {
    start_h: f64,
    start_dh: f64,
    env: &'a Environment,
}

impl Ray<'_> {
    pub fn from_h_ang(env: &Environment, h: f64, ang: f64) -> Ray {
        let dh = ang.tan();
        Ray {
            start_h: h,
            start_dh: dh,
            env,
        }
    }

    fn state_at_dist(&self, dist: f64) -> RayState {
        let tgt_x = dist.abs();

        let mut state = RayState {
            x: 0.0,
            h: self.start_h,
            dh: if dist >= 0.0 {
                self.start_dh
            } else {
                -self.start_dh
            },
        };

        let mut integrator = RK4Integrator::new(5.0);
        while state.x < tgt_x {
            integrator.propagate_in_place(
                &mut state,
                |state| self.env.calc_derivative_flat(state),
                StepSize::UseDefault,
            );
        }

        state
    }
}

impl<'a, 'b: 'a> Path<'a> for Ray<'b> {
    fn h_at_dist(&self, dist: f64) -> f64 {
        let state = self.state_at_dist(dist);
        state.h
    }

    fn angle_at_dist(&self, dist: f64) -> f64 {
        let state = self.state_at_dist(dist);
        state.get_angle(&self.env)
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
            |state| env.calc_derivative_flat(state),
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
