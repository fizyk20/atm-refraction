use path::Path;

#[derive(Clone, Copy, Debug)]
pub struct Line {
    pub a: f64,
    pub b: f64,
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
    fn start_h(&self) -> f64 {
        self.b
    }

    fn start_angle(&self) -> f64 {
        self.a.atan()
    }

    fn h_at_dist(&self, dist: f64) -> f64 {
        self.a * dist + self.b
    }

    fn angle_at_dist(&self, _dist: f64) -> f64 {
        self.a.atan()
    }
}
