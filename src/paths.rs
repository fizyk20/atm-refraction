pub trait Path {
    fn h_at_dist(&self, dist: f64) -> f64;
    fn angle_at_dist(&self, dist: f64) -> f64;
}

pub mod spherical {
    use super::Path;
    use crate::{Environment, RayState};
    use na::integration::{Integrator, RK4Integrator, StepSize};

    pub struct Line<'a> {
        env: &'a Environment,
        rmin: f64,
        phimin: f64,
    }

    impl Line<'_> {
        pub fn from_h_ang(env: &Environment, h: f64, ang: f64) -> Line {
            Line {
                env,
                rmin: (h + env.radius().unwrap()) * ang.cos(),
                phimin: -ang,
            }
        }

        pub fn from_two_points(env: &Environment, h1: f64, phi1: f64, h2: f64, phi2: f64) -> Line {
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

    impl Path for Line<'_> {
        fn h_at_dist(&self, dist: f64) -> f64 {
            let r = self.env.radius().unwrap();
            self.r(dist * 1e3 / r) - r
        }

        fn angle_at_dist(&self, dist: f64) -> f64 {
            dist * 1e3 / self.env.radius().unwrap() - self.phimin
        }
    }

    pub struct Ray<'a> {
        env: &'a Environment,
        start_h: f64,
        start_dh: f64,
    }

    impl Ray<'_> {
        pub fn from_h_ang(env: &Environment, h: f64, ang: f64) -> Ray {
            let dh = (h + env.radius().unwrap()) * ang.tan();
            Ray {
                env,
                start_h: h,
                start_dh: dh,
            }
        }

        fn state_at_dist(&self, dist: f64) -> RayState {
            let r = self.env.radius().unwrap();
            let tgt_phi = if dist >= 0.0 {
                dist * 1e3 / r
            } else {
                -dist * 1e3 / r
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

            let mut integrator = RK4Integrator::new(1e-6);
            while state.x < tgt_phi {
                integrator.propagate_in_place(
                    &mut state,
                    |state| self.env.calc_derivative_spherical(state),
                    StepSize::UseDefault,
                );
            }

            state
        }
    }

    impl Path for Ray<'_> {
        fn h_at_dist(&self, dist: f64) -> f64 {
            let state = self.state_at_dist(dist);
            state.h
        }

        fn angle_at_dist(&self, dist: f64) -> f64 {
            let state = self.state_at_dist(dist);
            (state.dr / (state.h + self.env.radius().unwrap())).atan()
        }
    }
}

pub mod flat {
    use super::Path;
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

    impl Path for Line {
        fn h_at_dist(&self, dist: f64) -> f64 {
            self.a * dist + self.b
        }

        fn angle_at_dist(&self, _dist: f64) -> f64 {
            self.a.atan()
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
                    |state| self.env.calc_derivative_flat(state),
                    StepSize::UseDefault,
                );
            }

            state
        }
    }

    impl Path for Ray<'_> {
        fn h_at_dist(&self, dist: f64) -> f64 {
            let state = self.state_at_dist(dist);
            state.h
        }

        fn angle_at_dist(&self, dist: f64) -> f64 {
            let state = self.state_at_dist(dist);
            state.dr.atan()
        }
    }
}
