pub trait Path {
    fn start_h(&self) -> f64;
    fn start_angle(&self) -> f64;
    fn h_at(&self, dist: f64) -> f64;
    fn angle_at(&self, dist: f64) -> f64;
}
