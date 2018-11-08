use flat;
use params::{EarthShape, Environment, Params, RayDir};
use spherical;
use PI;

pub trait Path {
    fn start_h(&self) -> f64;
    fn start_angle(&self) -> f64;
    fn h_at_dist(&self, dist: f64) -> f64;
    fn angle_at_dist(&self, dist: f64) -> f64;
}

fn ray_from_shape_h_ang(env: &Environment, h0: f64, angle: f64) -> Box<Path> {
    match env.shape {
        EarthShape::Spherical { radius } => Box::new(spherical::Ray::from_h_ang(
            env.atmosphere.clone(),
            radius,
            h0,
            angle,
        )),
        EarthShape::Flat => Box::new(flat::Ray::from_h_ang(env.atmosphere.clone(), h0, angle)),
    }
}

fn line_from_shape_h_ang(shape: &EarthShape, h0: f64, angle: f64) -> Box<Path> {
    match *shape {
        EarthShape::Spherical { radius } => {
            Box::new(spherical::Line::from_h_ang(radius, h0, angle))
        }
        EarthShape::Flat => Box::new(flat::Line::from_h_ang(h0, angle)),
    }
}

fn path_from_h_ang(env: &Environment, straight: bool, h0: f64, angle: f64) -> Box<Path> {
    if straight {
        line_from_shape_h_ang(&env.shape, h0, angle)
    } else {
        ray_from_shape_h_ang(env, h0, angle)
    }
}

fn find_angle_to_target(env: &Environment, h0: f64, tgt_h: f64, tgt_dist: f64) -> f64 {
    let (mut min_ang, mut max_ang) = (-1.5, 1.5);

    while max_ang - min_ang > 0.000001 {
        let cur_ang = 0.5 * (min_ang + max_ang);
        let ray = ray_from_shape_h_ang(env, h0, cur_ang);
        let h = ray.h_at_dist(tgt_dist);
        if h > tgt_h {
            max_ang = cur_ang;
        } else {
            min_ang = cur_ang;
        }
    }

    0.5 * (min_ang + max_ang)
}

fn path_from_h_tgt(
    env: &Environment,
    straight: bool,
    h0: f64,
    tgt_h: f64,
    tgt_dist: f64,
) -> Box<Path> {
    if straight {
        match env.shape {
            EarthShape::Spherical { radius } => Box::new(spherical::Line::from_two_points(
                radius,
                h0,
                0.0,
                tgt_h,
                tgt_dist * 1e3 / radius,
            )),
            EarthShape::Flat => {
                Box::new(flat::Line::from_two_points(h0, 0.0, tgt_h, tgt_dist * 1e3))
            }
        }
    } else {
        let ang = find_angle_to_target(env, h0, tgt_h, tgt_dist);
        ray_from_shape_h_ang(env, h0, ang)
    }
}

pub fn create_path(params: &Params) -> Box<Path> {
    match params.ray.dir {
        RayDir::Angle(ang) => {
            if params.verbose {
                println!("Starting angle: {} degrees from horizontal", ang);
            }
            path_from_h_ang(
                &params.env,
                params.straight,
                params.ray.start_h,
                ang * PI / 180.0,
            )
        }
        RayDir::Target { h, dist } => {
            if params.verbose {
                println!("Hits {} m ASL at a distance of {} km", h, dist);
            }
            path_from_h_tgt(&params.env, params.straight, params.ray.start_h, h, dist)
        }
        RayDir::Horizon => {
            if params.verbose {
                println!("Find angle to horizon");
            }
            path_from_h_ang(&params.env, params.straight, 0.0, 0.0)
        }
    }
}
