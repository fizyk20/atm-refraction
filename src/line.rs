#[derive(Clone, Copy, Debug)]
pub struct Line {
    pub rmin: f64,
    pub phimin: f64,
}

impl Line {
    pub fn from_r_dr(r: f64, phi: f64, dr: f64) -> Line {
        let dphi = (dr / r).atan();
        Line {
            rmin: r * dphi.cos(),
            phimin: phi - dphi,
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
