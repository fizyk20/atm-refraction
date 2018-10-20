//! US Standard Atmosphere model 1976
//! http://nebula.wsimg.com/ab321c1edd4fa69eaa94b5e8e769b113?AccessKeyId=AF1D67CEBF3A194F66A3&disposition=0&alloworigin=1

// gravitational acceleration
const G0: f64 = 9.80665;

// gas constant for air
const R: f64 = 287.053;

// lapse rates
const L0: f64 = -0.0065;
const L2: f64 = 0.001;
const L3: f64 = 0.0028;
const L5: f64 = -0.0028;
const L6: f64 = -0.002;

// base altitudes
const H1: f64 = 11e3;
const H2: f64 = 20e3;
const H3: f64 = 32e3;
const H4: f64 = 47e3;
const H5: f64 = 51e3;
const H6: f64 = 71e3;
const H7: f64 = 84.852e3;

#[allow(unused)]
pub fn density(p: f64, t: f64) -> f64 {
    p / (R * t)
}

// layer 0
#[inline]
fn temp0(t0: f64, h: f64) -> f64 {
    1.0 + h * L0 / t0
}

#[inline]
fn pressure0(t0: f64, h: f64) -> f64 {
    temp0(t0, h).powf(-G0 / L0 / R)
}

// layer 1
#[inline]
fn temp1(t0: f64, _h: f64) -> f64 {
    temp0(t0, H1)
}

#[inline]
fn pressure1(t0: f64, h: f64) -> f64 {
    let p1 = pressure0(t0, H1);
    let t1 = temp0(t0, H1) * t0;
    p1 * (-(h - H1) * G0 / R / t1).exp()
}

// layer 2
#[inline]
fn temp2(t0: f64, h: f64) -> f64 {
    temp1(t0, H2) + (h - H2) * L2 / t0
}

#[inline]
fn pressure2(t0: f64, h: f64) -> f64 {
    let t = temp2(t0, h) / temp1(t0, H2);
    pressure1(t0, H2) * t.powf(-G0 / L2 / R)
}

// layer 3
#[inline]
fn temp3(t0: f64, h: f64) -> f64 {
    temp2(t0, H3) + (h - H3) * L3 / t0
}

#[inline]
fn pressure3(t0: f64, h: f64) -> f64 {
    let t = temp3(t0, h) / temp2(t0, H3);
    pressure2(t0, H3) * t.powf(-G0 / L3 / R)
}

// layer 4
#[inline]
fn temp4(t0: f64, _h: f64) -> f64 {
    temp3(t0, H4)
}

#[inline]
fn pressure4(t0: f64, h: f64) -> f64 {
    let p4 = pressure3(t0, H4);
    let t4 = temp3(t0, H4) * t0;
    p4 * (-(h - H4) * G0 / R / t4).exp()
}

// layer 5
#[inline]
fn temp5(t0: f64, h: f64) -> f64 {
    temp4(t0, H5) + (h - H5) * L5 / t0
}

#[inline]
fn pressure5(t0: f64, h: f64) -> f64 {
    let t = temp5(t0, h) / temp4(t0, H5);
    pressure4(t0, H5) * t.powf(-G0 / L5 / R)
}

// layer 6
#[inline]
fn temp6(t0: f64, h: f64) -> f64 {
    temp5(t0, H6) + (h - H6) * L6 / t0
}

#[inline]
fn pressure6(t0: f64, h: f64) -> f64 {
    let t = temp6(t0, h) / temp5(t0, H6);
    pressure5(t0, H6) * t.powf(-G0 / L6 / R)
}

// layer 7
#[inline]
fn temp7(t0: f64, _h: f64) -> f64 {
    temp6(t0, H7)
}

#[inline]
fn pressure7(t0: f64, h: f64) -> f64 {
    let p7 = pressure6(t0, H7);
    let t7 = temp6(t0, H7) * t0;
    p7 * (-(h - H7) * G0 / R / t7).exp()
}

// Public functions
pub fn pressure(p0: f64, t0: f64, h: f64) -> f64 {
    if h < H1 {
        p0 * pressure0(t0, h)
    } else if h < H2 {
        p0 * pressure1(t0, h)
    } else if h < H3 {
        p0 * pressure2(t0, h)
    } else if h < H4 {
        p0 * pressure3(t0, h)
    } else if h < H5 {
        p0 * pressure4(t0, h)
    } else if h < H6 {
        p0 * pressure5(t0, h)
    } else if h < H7 {
        p0 * pressure6(t0, h)
    } else {
        p0 * pressure7(t0, h)
    }
}

pub fn temperature(t0: f64, h: f64) -> f64 {
    if h < H1 {
        t0 * temp0(t0, h)
    } else if h < H2 {
        t0 * temp1(t0, h)
    } else if h < H3 {
        t0 * temp2(t0, h)
    } else if h < H4 {
        t0 * temp3(t0, h)
    } else if h < H5 {
        t0 * temp4(t0, h)
    } else if h < H6 {
        t0 * temp5(t0, h)
    } else if h < H7 {
        t0 * temp6(t0, h)
    } else {
        t0 * temp7(t0, h)
    }
}
