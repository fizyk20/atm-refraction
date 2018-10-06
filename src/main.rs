extern crate numeric_algs as na;

mod ray;

use na::integration::{Integrator, RK4Integrator, StepSize};
use ray::{RayState, RayStateDerivative};

static R: f64 = 6_378_000.0;

#[inline]
fn n(h: f64) -> f64 {
    let n0 = 0.000293;
    let alpha = 1.25e-4;
    1.0 + n0 * (-alpha * h).exp()
}

#[inline]
fn dn(h: f64) -> f64 {
    let alpha = 1.25e-4;
    let n0 = 0.000293;
    -alpha * n0 * (-alpha * h).exp()
}

fn calc_derivative(state: &RayState) -> RayStateDerivative {
    let dr = state.dr;
    let h = state.h;

    let nr = n(h);
    let dnr = dn(h);

    let r = h + R;
    let d2r = dr * dr * dnr / nr + r * r * dnr / nr + 2.0 * dr * dr / r + r;

    RayStateDerivative {
        dphi: 1.0,
        dr: state.dr,
        d2r,
    }
}

#[derive(Clone, Copy, Debug)]
struct Line {
    rmin: f64,
    phimin: f64,
}

impl Line {
    fn from_r_dr(r: f64, phi: f64, dr: f64) -> Line {
        let dphi = (dr / r).atan();
        Line {
            rmin: r * dphi.cos(),
            phimin: phi - dphi,
        }
    }

    fn from_two_points(r1: f64, phi1: f64, r2: f64, phi2: f64) -> Line {
        let a = r1 / r2;
        let tanphi = (a * phi1.cos() - phi2.cos()) / (phi2.sin() - a * phi1.sin());
        let phimin = tanphi.atan();
        Line {
            rmin: r1 * (phi1 - phimin).cos(),
            phimin,
        }
    }

    fn r(&self, phi: f64) -> f64 {
        self.rmin / (phi - self.phimin).cos()
    }
}

fn find_height_at(dh0: f64, dist: f64) -> f64 {
    let h0 = 1565.0;

    let mut state = RayState {
        phi: 0.0,
        h: h0,
        dr: dh0,
    };

    let mut integrator = RK4Integrator::new(1e-5);
    while state.phi < dist / R {
        integrator.propagate_in_place(&mut state, calc_derivative, StepSize::UseDefault);
    }

    state.h
}

fn main() {
    let h0 = 1565.0;
    let (mut min_dh0, mut max_dh0) = (-280000.0, 0.0);

    while max_dh0 - min_dh0 > 0.001 {
        let cur_dh0 = 0.5 * (min_dh0 + max_dh0);
        let h = find_height_at(cur_dh0, 73e3);
        if h > 680.0 {
            max_dh0 = cur_dh0;
        } else {
            min_dh0 = cur_dh0;
        }
    }

    let dh0 = 0.5 * (min_dh0 + max_dh0);

    println!("Found dh0: {}", dh0);

    let h_at_73 = find_height_at(dh0, 73e3);
    let h_at_schneeberg = find_height_at(dh0, 277e3);

    println!("Altitude 73km from Praded: {}", h_at_73);
    println!("Altitude 277km from Praded: {}", h_at_schneeberg);
}
