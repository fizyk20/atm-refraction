use path::Path;

#[derive(Clone, Copy, Debug)]
pub struct Line {
    radius: f64,
    pub rmin: f64,
    pub phimin: f64,
}

impl Line {
    pub fn from_h_ang(radius: f64, h: f64, ang: f64) -> Line {
        Line {
            radius,
            rmin: (h + radius) * ang.cos(),
            phimin: -ang,
        }
    }

    pub fn from_two_points(radius: f64, r1: f64, phi1: f64, r2: f64, phi2: f64) -> Line {
        let a = r1 / r2;
        let tanphi = (a * phi1.cos() - phi2.cos()) / (phi2.sin() - a * phi1.sin());
        let phimin = tanphi.atan();
        Line {
            radius,
            rmin: r1 * (phi1 - phimin).cos(),
            phimin,
        }
    }

    pub fn r(&self, phi: f64) -> f64 {
        self.rmin / (phi - self.phimin).cos()
    }
}

impl Path for Line {
    fn start_h(&self) -> f64 {
        self.r(0.0) - self.radius
    }

    fn start_angle(&self) -> f64 {
        -self.phimin
    }

    fn h_at_dist(&self, dist: f64) -> f64 {
        self.r(dist * 1e3 / self.radius) - self.radius
    }

    fn angle_at_dist(&self, dist: f64) -> f64 {
        dist * 1e3 / self.radius - self.phimin
    }
}
