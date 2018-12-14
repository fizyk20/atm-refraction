use crate::Environment;
use na::{State, StateDerivative};
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug)]
pub struct RayState {
    pub x: f64,
    pub h: f64,
    pub dh: f64,
}

impl RayState {
    pub fn get_angle(&self, env: &Environment) -> f64 {
        if let Some(r) = env.radius() {
            (self.dh * r / (self.h + r)).atan()
        } else {
            self.dh.atan()
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RayStateDerivative {
    pub dx: f64,
    pub dh: f64,
    pub d2h: f64,
}

impl Add<RayStateDerivative> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn add(self, other: RayStateDerivative) -> RayStateDerivative {
        RayStateDerivative {
            dx: self.dx + other.dx,
            dh: self.dh + other.dh,
            d2h: self.d2h + other.d2h,
        }
    }
}

impl Sub<RayStateDerivative> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn sub(self, other: RayStateDerivative) -> RayStateDerivative {
        RayStateDerivative {
            dx: self.dx - other.dx,
            dh: self.dh - other.dh,
            d2h: self.d2h - other.d2h,
        }
    }
}

impl Mul<f64> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn mul(self, other: f64) -> RayStateDerivative {
        RayStateDerivative {
            dx: self.dx * other,
            dh: self.dh * other,
            d2h: self.d2h * other,
        }
    }
}

impl Div<f64> for RayStateDerivative {
    type Output = RayStateDerivative;
    fn div(self, other: f64) -> RayStateDerivative {
        RayStateDerivative {
            dx: self.dx / other,
            dh: self.dh / other,
            d2h: self.d2h / other,
        }
    }
}

impl Neg for RayStateDerivative {
    type Output = RayStateDerivative;
    fn neg(self) -> RayStateDerivative {
        RayStateDerivative {
            dx: -self.dx,
            dh: -self.dh,
            d2h: -self.d2h,
        }
    }
}

impl StateDerivative for RayStateDerivative {
    fn abs(&self) -> f64 {
        (self.dx * self.dx + self.dh * self.dh + self.d2h * self.d2h).sqrt()
    }
}

impl State for RayState {
    type Derivative = RayStateDerivative;
    fn shift_in_place(&mut self, dir: &RayStateDerivative, amount: f64) {
        self.x += dir.dx * amount;
        self.h += dir.dh * amount;
        self.dh += dir.d2h * amount;
    }
}
