use na::integration::{Integrator, RK4Integrator, StepSize};
use na::{State, StateDerivative};
use path::Path;
use std::ops::{Add, Div, Mul, Neg, Sub};

pub struct Ray {
    radius: f64,
    start_h: f64,
    start_dh: f64,
}

#[inline]
fn n(h: f64) -> f64 {
    let n0 = 0.000293;
    let alpha = 1.25e-4;
    1.0 + n0 * (-alpha * h).exp()
}

#[inline]
fn dn(h: f64) -> f64 {
    let alpha = 1.25e-4;
    let n0 = 0.000293;
    -alpha * n0 * (-alpha * h).exp()
}

fn calc_derivative(state: &RayState) -> RayStateDerivative {
    let dr = state.dr;
    let h = state.h;

    let nr = n(h);
    let dnr = dn(h);

    let r = h + state.R;
    let d2r = dr * dr * dnr / nr + r * r * dnr / nr + 2.0 * dr * dr / r + r;

    RayStateDerivative {
        dphi: 1.0,
        dr: state.dr,
        d2r,
    }
}

impl Ray {
    pub fn from_h_ang(radius: f64, h: f64, ang: f64) -> Ray {
        let dh = (h + radius) * ang.tan();
        Ray {
            radius,
            start_h: h,
            start_dh: dh,
        }
    }

    pub fn from_h_dh(radius: f64, h: f64, dh: f64) -> Ray {
        Ray {
            radius,
            start_h: h,
            start_dh: dh,
        }
    }

    fn state_at_phi(&self, phi: f64) -> RayState {
        let tgt_phi = if phi >= 0.0 { phi } else { -phi };
        let mut state = RayState {
            R: self.radius,
            phi: 0.0,
            h: self.start_h,
            dr: if phi >= 0.0 {
                self.start_dh
            } else {
                -self.start_dh
            },
        };

        let mut integrator = RK4Integrator::new(1e-5);
        while state.phi < tgt_phi {
            integrator.propagate_in_place(&mut state, calc_derivative, StepSize::UseDefault);
        }

        state
    }
}

impl Path for Ray {
    fn start_r(&self) -> f64 {
        self.start_h + self.radius
    }

    fn start_angle(&self) -> f64 {
        (self.start_dh / (self.start_h + self.radius)).atan()
    }

    fn r_at_phi(&self, phi: f64) -> f64 {
        let state = self.state_at_phi(phi);
        state.h + self.radius
    }

    fn angle_at_phi(&self, phi: f64) -> f64 {
        let state = self.state_at_phi(phi);
        (state.dr / (state.h + self.radius)).atan()
    }
}

#[derive(Clone, Copy, Debug)]
#[allow(non_snake_case)]
pub struct RayState {
    R: f64,
    pub phi: f64,
    pub h: f64,
    pub dr: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct RayStateDerivative {
    pub dphi: f64,
    pub dr: f64,
    pub d2r: f64,
}

impl Add<RayStateDerivative> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn add(self, other: RayStateDerivative) -> RayStateDerivative {
        RayStateDerivative {
            dphi: self.dphi + other.dphi,
            dr: self.dr + other.dr,
            d2r: self.d2r + other.d2r,
        }
    }
}

impl Sub<RayStateDerivative> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn sub(self, other: RayStateDerivative) -> RayStateDerivative {
        RayStateDerivative {
            dphi: self.dphi - other.dphi,
            dr: self.dr - other.dr,
            d2r: self.d2r - other.d2r,
        }
    }
}

impl Mul<f64> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn mul(self, other: f64) -> RayStateDerivative {
        RayStateDerivative {
            dphi: self.dphi * other,
            dr: self.dr * other,
            d2r: self.d2r * other,
        }
    }
}

impl Div<f64> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn div(self, other: f64) -> RayStateDerivative {
        RayStateDerivative {
            dphi: self.dphi / other,
            dr: self.dr / other,
            d2r: self.d2r / other,
        }
    }
}

impl Neg for RayStateDerivative {
    type Output = RayStateDerivative;
    fn neg(self) -> RayStateDerivative {
        RayStateDerivative {
            dphi: -self.dphi,
            dr: -self.dr,
            d2r: -self.d2r,
        }
    }
}

impl StateDerivative for RayStateDerivative {
    fn abs(&self) -> f64 {
        (self.dphi * self.dphi + self.dr * self.dr + self.d2r * self.d2r).sqrt()
    }
}

impl State for RayState {
    type Derivative = RayStateDerivative;
    fn shift_in_place(&mut self, dir: &RayStateDerivative, amount: f64) {
        self.phi += dir.dphi * amount;
        self.h += dir.dr * amount;
        self.dr += dir.d2r * amount;
    }
}
