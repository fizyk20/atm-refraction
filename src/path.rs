pub trait Path {
    fn start_r(&self) -> f64;
    fn start_angle(&self) -> f64;
    fn r_at_phi(&self, phi: f64) -> f64;
    fn angle_at_phi(&self, phi: f64) -> f64;
}
