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

pub static PI: f64 = 3.1415926536;

fn find_ray_from_target(radius: f64, h0: f64, tgt_h: f64, tgt_dist: f64) -> Ray {
    let (mut min_ang, mut max_ang) = (-1.5, 1.5);

    while max_ang - min_ang > 0.000001 {
        let cur_ang = 0.5 * (min_ang + max_ang);
        let ray = Ray::from_h_ang(radius, h0, cur_ang);
        let h = ray.r_at_phi(tgt_dist * 1e3 / radius) - radius;
        if h > tgt_h {
            max_ang = cur_ang;
        } else {
            min_ang = cur_ang;
        }
    }

    Ray::from_h_ang(radius, h0, 0.5 * (min_ang + max_ang))
}

fn find_dist_for_h(radius: f64, ray: &Path, tgt_h: f64) -> f64 {
    let (mut min_dist, mut max_dist) = (0.0, 5000.0);

    while max_dist - min_dist > 0.00001 {
        let cur_dist = 0.5 * (min_dist + max_dist);
        let h = ray.r_at_phi(cur_dist * 1e3 / radius) - radius;
        if h > tgt_h {
            max_dist = cur_dist;
        } else {
            min_dist = cur_dist;
        }
    }

    0.5 * (min_dist + max_dist)
}

fn main() {
    let params = parse_arguments();

    if params.verbose {
        println!("Ray parameters chosen:");
        println!("Earth radius: {} km", params.radius / 1e3);
        println!("Starting altitude: {} m ASL", params.ray.start_h);
    }

    let ray: Box<Path> = match params.ray.dir {
        RayDir::Angle(ang) => {
            if params.verbose {
                println!("Starting angle: {} degrees from horizontal", ang);
            }
            if params.straight {
                Box::new(Line::from_r_ang(
                    params.ray.start_h + params.radius,
                    ang * PI / 180.0,
                ))
            } else {
                Box::new(Ray::from_h_ang(
                    params.radius,
                    params.ray.start_h,
                    ang * PI / 180.0,
                ))
            }
        }
        RayDir::Target { h, dist } => {
            if params.verbose {
                println!("Hits {} m ASL at a distance of {} km", h, dist);
            }
            if params.straight {
                Box::new(Line::from_two_points(
                    params.ray.start_h + params.radius,
                    0.0,
                    h + params.radius,
                    dist * 1e3 / params.radius,
                ))
            } else {
                Box::new(find_ray_from_target(
                    params.radius,
                    params.ray.start_h,
                    h,
                    dist,
                ))
            }
        }
        RayDir::Horizon => {
            if params.verbose {
                println!("Find angle to horizon");
            }
            if params.straight {
                Box::new(Line::from_r_ang(params.radius, 0.0))
            } else {
                Box::new(Ray::from_h_dh(params.radius, 0.0, 0.0))
            }
        }
    };
    if params.straight && params.verbose {
        println!("Straight-line calculation chosen.");
    }
    if params.verbose {
        println!();
    }

    for output in &params.output {
        match *output {
            Output::HAtDist(dist) => {
                if params.verbose {
                    println!(
                        "Altitude at distance {} km: {}",
                        dist,
                        ray.r_at_phi(dist * 1e3 / params.radius) - params.radius
                    );
                } else {
                    println!(
                        "{}",
                        ray.r_at_phi(dist * 1e3 / params.radius) - params.radius
                    );
                }
            }
            Output::Angle => {
                if params.verbose {
                    println!("Starting angle: {} degrees", ray.start_angle() * 180.0 / PI);
                } else {
                    println!("{}", ray.start_angle() * 180.0 / PI);
                }
            }
            Output::Horizon => {
                let dist_to_target_h = find_dist_for_h(params.radius, &*ray, params.ray.start_h);
                let ang = ray.angle_at_phi(dist_to_target_h * 1e3 / params.radius);
                if params.verbose {
                    println!("Angle to the horizon: {} degrees", -ang * 180.0 / PI);
                } else {
                    println!("{}", -ang * 180.0 / PI);
                }
            }
        }
    }
}
