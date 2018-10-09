use clap::{App, Arg};

/// Ray direction description
pub enum RayDir {
    /// angle from the horizon
    Angle(f64),
    /// hit a given altitude at the given distance
    Target { h: f64, dist: f64 },
    /// special value for finding the horizon
    Horizon,
}

/// Definition of a ray
pub struct RayData {
    /// starting altitude
    pub start_h: f64,
    /// direction of propagation
    pub dir: RayDir,
}

/// what info to output
pub enum Output {
    /// altitude at a given distance
    HAtDist(f64),
    /// output starting angle of the ray
    Angle,
    /// output the angle to the horizon
    Horizon,
}

pub struct Params {
    pub ray: RayData,
    pub radius: f64,
    pub straight: bool,
    pub output: Vec<Output>,
}

pub fn parse_arguments() -> Params {
    let matches = App::new("Atmospheric Refraction Calculator")
        .about("Calculates paths of light in Earth's atmosphere")
        .arg(
            Arg::with_name("start_h")
                .short("h")
                .long("start-h")
                .value_name("ALTITUDE")
                .help("Starting point altitude (meters) (default = 1)")
                .takes_value(true),
        ).arg(
            Arg::with_name("start_angle")
                .short("a")
                .long("start-angle")
                .value_name("ANGLE")
                .help("Starting direction, angle relative to horizontal (degrees); conflicts with --tgt-h and --tgt-dist")
                .takes_value(true),
        ).arg(
            Arg::with_name("target_h")
                .short("t")
                .long("tgt-h")
                .value_name("ALTITUDE")
                .help("Target point altitude (meters); conflicts with --start-angle")
                .takes_value(true),
        ).arg(
            Arg::with_name("target_dist")
                .short("d")
                .long("tgt-dist")
                .value_name("DISTANCE")
                .help("Target point distance (kilometers); conflicts with --start-angle")
                .takes_value(true),
        ).arg(
            Arg::with_name("radius")
                .short("R")
                .long("radius")
                .value_name("RADIUS")
                .help("Calculate assuming the given value as the Earth's radius, in km (default: 6378)")
                .takes_value(true),
        ).arg(
            Arg::with_name("output_dist")
                .short("o")
                .long("output-dist")
                .value_name("DISTANCE")
                .help("Distance at which to output altitude (kilometers)")
                .takes_value(true),
        ).arg(
            Arg::with_name("output_ang")
                .long("output-ang")
                .help("Output the starting angle of the ray")
                .takes_value(false),
        ).arg(
            Arg::with_name("output_horizon")
                .long("output-horizon")
                .help("Output the angle to the horizon")
                .takes_value(false),
        ).arg(
            Arg::with_name("straight")
                .short("s")
                .long("straight")
                .help("Calculation for a straight-line ray")
                .takes_value(false),
        ).get_matches();
    let start_h: f64 = matches
        .value_of("start_h")
        .unwrap_or("1.0")
        .parse()
        .ok()
        .expect("Invalid altitude passed to start-h");
    let radius: f64 = matches
        .value_of("radius")
        .unwrap_or("6378.0")
        .parse()
        .ok()
        .expect("Invalid radius passed");
    let start_angle = matches.value_of("start_angle");
    let tgt_h = matches.value_of("target_h");
    let tgt_dist = matches.value_of("target_dist");

    let ray_dir = if matches.is_present("output_horizon") {
        RayDir::Horizon
    } else {
        match (start_angle, tgt_h, tgt_dist) {
            (Some(ang), None, None) => RayDir::Angle(
                ang.parse()
                    .ok()
                    .expect("Invalid angle passed to --start-angle"),
            ),
            (None, Some(h), Some(dist)) => RayDir::Target {
                h: h.parse().ok().expect("Invalid altitude passed to --tgt-h"),
                dist: dist
                    .parse()
                    .ok()
                    .expect("Invalid distance passed to --tgt-dist"),
            },
            (None, None, None) => panic!("No ray direction chosen!"),
            _ => panic!("Conflicting options detected (--start-angle, --tgt-h, --tgt-dist)"),
        }
    };
    let ray = RayData {
        start_h,
        dir: ray_dir,
    };
    let mut output = Vec::new();
    if let Some(dist) = matches
        .value_of("output_dist")
        .and_then(|val| val.parse().ok())
    {
        output.push(Output::HAtDist(dist));
    }
    if matches.is_present("output_ang") {
        output.push(Output::Angle);
    }
    if matches.is_present("output_horizon") {
        output = vec![Output::Horizon];
    }
    Params {
        ray,
        straight: matches.is_present("straight"),
        radius: radius * 1e3,
        output,
    }
}
