use air::{air_index_minus_1, pressure, temperature};
use na::{State, StateDerivative};
use std::ops::{Add, Div, Mul, Neg, Sub};

#[inline]
fn n_minus_1(h: f64) -> f64 {
    let p0 = 101325.0;
    let t0 = 273.15;

    let pressure = pressure(p0, t0, h);
    let temperature = temperature(t0, h);
    let rh = 0.0;

    air_index_minus_1(530e-9, pressure, temperature, rh)
}

#[inline]
pub fn n(h: f64) -> f64 {
    n_minus_1(h) + 1.0
}

#[inline]
fn dn(h: f64) -> f64 {
    let epsilon = 0.01;
    let n1 = n_minus_1(h - epsilon);
    let n2 = n_minus_1(h + epsilon);
    (n2 - n1) / (2.0 * epsilon)
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
