use path::Path;

#[derive(Clone, Copy, Debug)]
pub struct Line {
    pub rmin: f64,
    pub phimin: f64,
}

impl Line {
    pub fn from_r_ang(r: f64, ang: f64) -> Line {
        Line {
            rmin: r * ang.cos(),
            phimin: -ang,
        }
    }

    pub fn from_two_points(r1: f64, phi1: f64, r2: f64, phi2: f64) -> Line {
        let a = r1 / r2;
        let tanphi = (a * phi1.cos() - phi2.cos()) / (phi2.sin() - a * phi1.sin());
        let phimin = tanphi.atan();
        Line {
            rmin: r1 * (phi1 - phimin).cos(),
            phimin,
        }
    }

    pub fn r(&self, phi: f64) -> f64 {
        self.rmin / (phi - self.phimin).cos()
    }
}

impl Path for Line {
    fn start_r(&self) -> f64 {
        self.r(0.0)
    }

    fn start_angle(&self) -> f64 {
        -self.phimin
    }

    fn r_at_phi(&self, phi: f64) -> f64 {
        self.r(phi)
    }

    fn angle_at_phi(&self, phi: f64) -> f64 {
        phi - self.phimin
    }
}
