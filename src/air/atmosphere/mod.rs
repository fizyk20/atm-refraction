mod pressure_profile;
pub mod vertical_profile;

use self::{
    pressure_profile::PressureProfile,
    vertical_profile::{FunctionDef, VerticalProfile, VerticalProfileBuilder},
};

#[cfg(feature = "serialization")]
use cubic_splines::BoundaryCondition;

/// mu*g/R
pub const A: f64 = 0.03416320331088684;

#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct PressureFixedPoint {
    altitude: f64,
    pressure: f64,
}

#[derive(Clone, Debug)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct FunctionDefWithAlt {
    altitude: f64,
    function: FunctionDef,
}

#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct TemperatureFixedPoint {
    altitude: f64,
    temperature: f64,
}

#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct HumidityFixedPoint {
    altitude: f64,
    humidity: f64,
}

#[derive(Clone, Debug)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct AtmosphereDef {
    #[cfg_attr(feature = "serialization", serde(default = "default_pressure"))]
    pressure: PressureFixedPoint,
    first_temperature_function: FunctionDef,
    #[cfg_attr(feature = "serialization", serde(default))]
    next_functions: Vec<FunctionDefWithAlt>,
    temperature_fixed_point: Option<TemperatureFixedPoint>,

    #[cfg_attr(
        feature = "serialization",
        serde(default = "default_first_humidity_function")
    )]
    first_humidity_function: FunctionDef,
    #[cfg_attr(feature = "serialization", serde(default))]
    next_humidity_functions: Vec<FunctionDefWithAlt>,
    humidity_fixed_point: Option<HumidityFixedPoint>,
}

impl AtmosphereDef {
    pub fn us_76() -> Self {
        AtmosphereDef {
            pressure: PressureFixedPoint {
                altitude: 0.0,
                pressure: 101325.0,
            },
            first_temperature_function: FunctionDef::Linear { gradient: -0.0065 },
            next_functions: vec![
                FunctionDefWithAlt {
                    altitude: 11e3,
                    function: FunctionDef::Linear { gradient: 0.0 },
                },
                FunctionDefWithAlt {
                    altitude: 20e3,
                    function: FunctionDef::Linear { gradient: 0.001 },
                },
                FunctionDefWithAlt {
                    altitude: 32e3,
                    function: FunctionDef::Linear { gradient: 0.0028 },
                },
                FunctionDefWithAlt {
                    altitude: 47e3,
                    function: FunctionDef::Linear { gradient: 0.0 },
                },
                FunctionDefWithAlt {
                    altitude: 51e3,
                    function: FunctionDef::Linear { gradient: -0.0028 },
                },
                FunctionDefWithAlt {
                    altitude: 71e3,
                    function: FunctionDef::Linear { gradient: -0.002 },
                },
                FunctionDefWithAlt {
                    altitude: 84.852e3,
                    function: FunctionDef::Linear { gradient: 0.0 },
                },
            ],
            temperature_fixed_point: Some(TemperatureFixedPoint {
                altitude: 0.0,
                temperature: 288.0,
            }),
            first_humidity_function: FunctionDef::Linear { gradient: 0.0 },
            next_humidity_functions: vec![],
            humidity_fixed_point: Some(HumidityFixedPoint {
                altitude: 0.0,
                humidity: 0.0,
            }),
        }
    }
}

#[cfg(feature = "serialization")]
fn default_pressure() -> PressureFixedPoint {
    PressureFixedPoint {
        altitude: 0.0,
        pressure: 101325.0,
    }
}

#[cfg(feature = "serialization")]
fn default_first_humidity_function() -> FunctionDef {
    FunctionDef::Spline {
        points: vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)],
        boundary_condition: BoundaryCondition::Natural,
    }
}

/// A structure representing an atmospheric model. It provides the temperature and density as
/// functions of altitude
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct Atmosphere {
    pressure: PressureProfile,
    temperature: VerticalProfile,
    humidity: VerticalProfile,
}

impl Atmosphere {
    /// Creates the atmospheric model from a parsed definition.
    pub fn from_def(def: AtmosphereDef) -> Atmosphere {
        let mut builder = VerticalProfileBuilder::new(def.first_temperature_function);
        if let Some(point) = def.temperature_fixed_point {
            builder = builder.with_fixed_value(point.altitude, point.temperature);
        }
        for fun_def in def.next_functions {
            builder = builder.with_next_function(fun_def.altitude, fun_def.function);
        }
        let temperature = builder.build().unwrap();

        let mut builder = VerticalProfileBuilder::new(def.first_humidity_function);
        if let Some(point) = def.humidity_fixed_point {
            builder = builder.with_fixed_value(point.altitude, point.humidity);
        }
        for fun_def in def.next_humidity_functions {
            builder = builder.with_next_function(fun_def.altitude, fun_def.function);
        }
        let humidity = builder.build().unwrap();

        let pressure = PressureProfile::from_temperature_profile(
            &temperature,
            def.pressure.pressure,
            def.pressure.altitude,
        );

        Atmosphere {
            pressure,
            temperature,
            humidity,
        }
    }

    /// Returns the temperature at the given altitude
    pub fn temperature(&self, h: f64) -> f64 {
        self.temperature.eval(h)
    }

    /// Returns the derivative of temperature with respect to altitude at the given altitude
    pub fn dtemperature(&self, h: f64) -> f64 {
        self.temperature.eval_derivative(h)
    }

    /// Returns the pressure at the given altitude
    pub fn pressure(&self, h: f64) -> f64 {
        self.pressure.eval(h)
    }

    /// Returns the derivative of pressure at the given altitude
    pub fn dpressure(&self, h: f64) -> f64 {
        let p = self.pressure(h);
        let t = self.temperature(h);
        -A * p / t
    }

    /// Returns the temperature at the given altitude
    pub fn humidity(&self, h: f64) -> f64 {
        self.humidity.eval(h)
    }

    /// Returns the derivative of temperature with respect to altitude at the given altitude
    pub fn dhumidity(&self, h: f64) -> f64 {
        self.humidity.eval_derivative(h)
    }
}

/// Returns the US-1976 standard model of the Earth's atmosphere.
///
/// The temperatures are expressed in kelvins (K), and the pressure in hectopascals (hPa).
pub fn us76_atmosphere() -> Atmosphere {
    let atm_def = AtmosphereDef::us_76();
    Atmosphere::from_def(atm_def)
}

#[cfg(test)]
mod test {
    use super::*;

    use cubic_splines::BoundaryCondition;

    #[test]
    fn test_us76() {
        let atmosphere = Atmosphere::from_def(AtmosphereDef::us_76());
        assert_eq!(atmosphere.pressure(0.0), 101325.0);
        assert_eq!(atmosphere.temperature(0.0), 288.0);
    }

    #[test]
    fn test_spline() {
        let atmosphere_def = AtmosphereDef {
            pressure: PressureFixedPoint {
                altitude: 0.0,
                pressure: 1000.0,
            },
            first_temperature_function: FunctionDef::Spline {
                boundary_condition: BoundaryCondition::Derivatives(-0.0065, -0.0065),
                points: vec![
                    (0.0, 281.6),
                    (12.5, 283.4),
                    (19.4, 281.9),
                    (24.0, 284.7),
                    (34.0, 290.5),
                ],
            },
            next_functions: vec![],
            temperature_fixed_point: None,
            ..AtmosphereDef::us_76()
        };
        let atmosphere = Atmosphere::from_def(atmosphere_def);
        for i in 0..600 {
            let h = i as f64 * 0.5;
            println!(
                "{} {} {}",
                h,
                atmosphere.temperature(h),
                atmosphere.pressure(h)
            );
        }
    }
}
