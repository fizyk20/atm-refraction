use na::{State, StateDerivative};
use std::ops::{Add, Div, Mul, Neg, Sub};

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

pub fn calc_derivative_spherical(radius: f64, state: &RayState) -> RayStateDerivative {
    let dr = state.dr;
    let h = state.h;

    let nr = n(h);
    let dnr = dn(h);

    let r = h + radius;
    let d2r = dr * dr * dnr / nr + r * r * dnr / nr + 2.0 * dr * dr / r + r;

    RayStateDerivative { dx: 1.0, dr, d2r }
}

pub fn calc_derivative_flat(state: &RayState) -> RayStateDerivative {
    let dr = state.dr;
    let h = state.h;

    let nr = n(h);
    let dnr = dn(h);

    let d2r = dnr / nr * (1.0 + dr * dr);

    RayStateDerivative { dx: 1.0, dr, d2r }
}

#[derive(Clone, Copy, Debug)]
pub struct RayState {
    pub x: f64,
    pub h: f64,
    pub dr: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct RayStateDerivative {
    pub dx: f64,
    pub dr: f64,
    pub d2r: f64,
}

impl Add<RayStateDerivative> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn add(self, other: RayStateDerivative) -> RayStateDerivative {
        RayStateDerivative {
            dx: self.dx + other.dx,
            dr: self.dr + other.dr,
            d2r: self.d2r + other.d2r,
        }
    }
}

impl Sub<RayStateDerivative> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn sub(self, other: RayStateDerivative) -> RayStateDerivative {
        RayStateDerivative {
            dx: self.dx - other.dx,
            dr: self.dr - other.dr,
            d2r: self.d2r - other.d2r,
        }
    }
}

impl Mul<f64> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn mul(self, other: f64) -> RayStateDerivative {
        RayStateDerivative {
            dx: self.dx * other,
            dr: self.dr * other,
            d2r: self.d2r * other,
        }
    }
}

impl Div<f64> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn div(self, other: f64) -> RayStateDerivative {
        RayStateDerivative {
            dx: self.dx / other,
            dr: self.dr / other,
            d2r: self.d2r / other,
        }
    }
}

impl Neg for RayStateDerivative {
    type Output = RayStateDerivative;
    fn neg(self) -> RayStateDerivative {
        RayStateDerivative {
            dx: -self.dx,
            dr: -self.dr,
            d2r: -self.d2r,
        }
    }
}

impl StateDerivative for RayStateDerivative {
    fn abs(&self) -> f64 {
        (self.dx * self.dx + self.dr * self.dr + self.d2r * self.d2r).sqrt()
    }
}

impl State for RayState {
    type Derivative = RayStateDerivative;
    fn shift_in_place(&mut self, dir: &RayStateDerivative, amount: f64) {
        self.x += dir.dx * amount;
        self.h += dir.dr * amount;
        self.dr += dir.d2r * amount;
    }
}
