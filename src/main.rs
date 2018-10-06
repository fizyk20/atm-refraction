extern crate clap;
extern crate numeric_algs as na;

mod line;
mod params;
mod ray;

use line::Line;
use na::integration::{Integrator, RK4Integrator, StepSize};
use params::*;
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

fn find_height_at(h0: f64, dh0: f64, dist: f64) -> f64 {
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

fn find_dr_from_target(h0: f64, tgt_h: f64, tgt_dist: f64) -> f64 {
    let (mut min_dh0, mut max_dh0) = (-1e9, 0.0);

    while max_dh0 - min_dh0 > 0.001 {
        let cur_dh0 = 0.5 * (min_dh0 + max_dh0);
        let h = find_height_at(h0, cur_dh0, tgt_dist * 1e3);
        if h > tgt_h {
            max_dh0 = cur_dh0;
        } else {
            min_dh0 = cur_dh0;
        }
    }

    0.5 * (min_dh0 + max_dh0)
}

fn main() {
    let params = parse_arguments();

    println!("Ray parameters chosen:");
    println!("Starting altitude: {} m ASL", params.ray.start_h);

    let dh0 = match params.ray.dir {
        RayDir::Angle(ang) => {
            println!("Starting angle: {} degrees from horizontal", ang);
            (params.ray.start_h + R) * (ang * 3.1415926535 / 180.0).tan()
        }
        RayDir::Target { h, dist } => {
            println!("Hits {} m ASL at a distance of {} km", h, dist);
            if params.straight {
                let line =
                    Line::from_two_points(params.ray.start_h + R, 0.0, h + R, dist * 1e3 / R);
                line.r(0.0) * (-line.phimin).tan()
            } else {
                find_dr_from_target(params.ray.start_h, h, dist)
            }
        }
    };
    if params.straight {
        println!("Straight-line calculation chosen.");
    }
    println!();

    match params.output {
        Output::HAtDist(dist) => {
            if params.straight {
                let line = Line::from_r_dr(params.ray.start_h + R, 0.0, dh0);
                println!(
                    "Straight-line altitude at distance {} km: {}",
                    dist,
                    line.r(dist * 1e3 / R) - R
                );
            } else {
                println!(
                    "Altitude at distance {} km: {}",
                    dist,
                    find_height_at(params.ray.start_h, dh0, dist * 1e3)
                );
            }
        }
    }
}
