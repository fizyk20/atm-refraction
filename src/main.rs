extern crate clap;
extern crate numeric_algs as na;

mod line;
mod params;
mod path;
mod ray;

use line::Line;
use params::*;
use path::Path;
use ray::Ray;

pub static R: f64 = 6_378_000.0;
pub static PI: f64 = 3.1415926536;

fn find_ray_from_target(h0: f64, tgt_h: f64, tgt_dist: f64) -> Ray {
    let (mut min_ang, mut max_ang) = (-1.5, 1.5);

    while max_ang - min_ang > 0.000001 {
        let cur_ang = 0.5 * (min_ang + max_ang);
        let ray = Ray::from_h_ang(h0, cur_ang);
        let h = ray.h_at(tgt_dist);
        if h > tgt_h {
            max_ang = cur_ang;
        } else {
            min_ang = cur_ang;
        }
    }

    Ray::from_h_ang(h0, 0.5 * (min_ang + max_ang))
}

fn main() {
    let params = parse_arguments();

    println!("Ray parameters chosen:");
    println!("Starting altitude: {} m ASL", params.ray.start_h);

    let ray: Box<Path> = match params.ray.dir {
        RayDir::Angle(ang) => {
            println!("Starting angle: {} degrees from horizontal", ang);
            if params.straight {
                Box::new(Line::from_r_ang(params.ray.start_h + R, ang * PI / 180.0))
            } else {
                Box::new(Ray::from_h_ang(params.ray.start_h, ang * PI / 180.0))
            }
        }
        RayDir::Target { h, dist } => {
            println!("Hits {} m ASL at a distance of {} km", h, dist);
            if params.straight {
                Box::new(Line::from_two_points(
                    params.ray.start_h + R,
                    0.0,
                    h + R,
                    dist * 1e3 / R,
                ))
            } else {
                Box::new(find_ray_from_target(params.ray.start_h, h, dist))
            }
        }
    };
    if params.straight {
        println!("Straight-line calculation chosen.");
    }
    println!();

    for output in &params.output {
        match *output {
            Output::HAtDist(dist) => {
                println!("Altitude at distance {} km: {}", dist, ray.h_at(dist));
            }
            Output::Angle => {
                println!("Starting angle: {} degrees", ray.start_angle() * 180.0 / PI);
            }
        }
    }
}
