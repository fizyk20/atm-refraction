use na::{State, StateDerivative};
use std::ops::{Add, Div, Mul, Neg, Sub};

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
