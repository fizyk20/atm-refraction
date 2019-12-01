use super::parser::{parse_atmosphere, AtmosphereDef, LapseDef, PressureDef, TemperatureAt};
use std::cmp::Ordering;
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// A structure representing an atmospheric model. It provides the temperature and density as
/// functions of altitude
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct Atmosphere {
    layer_altitudes: Vec<f64>,
    first_lapse: f64,
    lapses: Vec<f64>,
    start_pressures: Vec<f64>,
    start_temperatures: Vec<f64>,
}

fn cmp_f64(a: &f64, b: &f64) -> Ordering {
    if a < b {
        Ordering::Less
    } else if a > b {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

// mu*g/R
const A: f64 = 0.03416320331088684;

fn isothermal(t0: f64, dh: f64) -> f64 {
    (-A * dh / t0).exp()
}

fn non_isothermal(t0: f64, lapse: f64, dh: f64) -> f64 {
    (1.0 + lapse / t0 * dh).powf(-A / lapse)
}

fn shift_p_with_lapse(t0: f64, lapse: f64, dh: f64) -> f64 {
    if lapse == 0.0 {
        isothermal(t0, dh)
    } else {
        non_isothermal(t0, lapse, dh)
    }
}

impl Atmosphere {
    fn get_lapses(lapses: &[LapseDef]) -> (f64, Vec<f64>, Vec<f64>) {
        let mut first_lapses = lapses
            .into_iter()
            .filter(|lapse_def| lapse_def.start_h.is_none());
        let first_lapse = *first_lapses.next().unwrap();
        assert!(
            first_lapses.next().is_none(),
            "There was more than one lapse defined with no altitude!"
        );

        let mut lapses = lapses
            .iter()
            .filter_map(|lapse_def| lapse_def.start_h.map(|h| (h, lapse_def.lapse)))
            .collect::<Vec<_>>();
        lapses.sort_by(|a, b| cmp_f64(&a.0, &b.0));

        let altitudes = lapses.iter().map(|l| l.0).collect();
        let lapses = lapses.iter().map(|l| l.1).collect();

        (first_lapse.lapse, lapses, altitudes)
    }

    fn find_t0(h0: f64, th0: f64, first_lapse: f64, lapses: &[f64], alts: &[f64]) -> f64 {
        match alts.binary_search_by(|a| cmp_f64(a, &h0)) {
            Ok(mut i) => {
                let mut t0 = th0;
                while i > 0 {
                    t0 -= lapses[i - 1] * (alts[i] - alts[i - 1]);
                    i -= 1;
                }
                t0
            }
            Err(0) => th0 + first_lapse * (alts[0] - h0),
            Err(mut i) => {
                i -= 1;
                let mut t0 = th0 - lapses[i] * (h0 - alts[i]);
                while i > 0 {
                    t0 -= lapses[i - 1] * (alts[i] - alts[i - 1]);
                    i -= 1;
                }
                t0
            }
        }
    }

    fn find_p0(
        h0: f64,
        ph0: f64,
        first_lapse: f64,
        lapses: &[f64],
        alts: &[f64],
        temps: &[f64],
    ) -> f64 {
        match alts.binary_search_by(|a| cmp_f64(a, &h0)) {
            Ok(mut i) => {
                let mut p0 = ph0;
                while i > 0 {
                    p0 /= shift_p_with_lapse(temps[i - 1], lapses[i - 1], alts[i] - alts[i - 1]);
                    i -= 1;
                }
                p0
            }
            Err(0) => ph0 / shift_p_with_lapse(temps[0], first_lapse, h0 - alts[0]),
            Err(mut i) => {
                i -= 1;
                let mut p0 = ph0 / shift_p_with_lapse(temps[i], lapses[i], h0 - alts[i]);
                while i > 0 {
                    p0 /= shift_p_with_lapse(temps[i - 1], lapses[i - 1], alts[i] - alts[i - 1]);
                    i -= 1;
                }
                p0
            }
        }
    }

    /// Creates the atmospheric model from a parsed definition.
    pub fn from_def(def: AtmosphereDef) -> Atmosphere {
        let PressureDef {
            start_h: h0_p,
            start_p: p0,
        } = def.pressure;
        let (first_lapse, lapses, layer_altitudes) = Self::get_lapses(&def.temperature.lapses);
        let TemperatureAt {
            start_h: h0_t,
            start_t: t0,
        } = def.temperature.start;

        let t0 = Self::find_t0(h0_t, t0, first_lapse, &lapses, &layer_altitudes);
        let mut start_temperatures = vec![t0];
        let mut t = t0;
        for (i, alt) in layer_altitudes.iter().enumerate().skip(1) {
            t += lapses[i - 1] * (alt - layer_altitudes[i - 1]);
            start_temperatures.push(t);
        }

        let p0 = Self::find_p0(
            h0_p,
            p0,
            first_lapse,
            &lapses,
            &layer_altitudes,
            &start_temperatures,
        );
        let mut start_pressures = vec![p0];
        let mut p = p0;
        for (i, alt) in layer_altitudes.iter().enumerate().skip(1) {
            p *= shift_p_with_lapse(
                start_temperatures[i - 1],
                lapses[i - 1],
                alt - layer_altitudes[i - 1],
            );
            start_pressures.push(p);
        }

        Atmosphere {
            first_lapse,
            lapses,
            layer_altitudes,
            start_pressures,
            start_temperatures,
        }
    }

    /// Returns the temperature at the given altitude
    pub fn temperature(&self, h: f64) -> f64 {
        match self.layer_altitudes.binary_search_by(|a| cmp_f64(a, &h)) {
            Ok(i) => self.start_temperatures[i],
            Err(0) => self.start_temperatures[0] - self.first_lapse * (self.layer_altitudes[0] - h),
            Err(i) => {
                self.start_temperatures[i - 1]
                    + self.lapses[i - 1] * (h - self.layer_altitudes[i - 1])
            }
        }
    }

    /// Returns the pressure at the given altitude
    pub fn pressure(&self, h: f64) -> f64 {
        match self.layer_altitudes.binary_search_by(|a| cmp_f64(a, &h)) {
            Ok(i) => self.start_pressures[i],
            Err(0) => {
                self.start_pressures[0]
                    * shift_p_with_lapse(
                        self.start_temperatures[0],
                        self.first_lapse,
                        h - self.layer_altitudes[0],
                    )
            }
            Err(i) => {
                self.start_pressures[i - 1]
                    * shift_p_with_lapse(
                        self.start_temperatures[i - 1],
                        self.lapses[i - 1],
                        h - self.layer_altitudes[i - 1],
                    )
            }
        }
    }
}

/// Parses an atmosphere definition from a string
pub fn atm_from_str<'a>(def: &'a str) -> Result<Atmosphere, nom::Err<nom::types::CompleteStr<'a>>> {
    let atm_def = parse_atmosphere(def)?;
    Ok(Atmosphere::from_def(atm_def.1))
}

/// Reads an atmosphere definition from file and returns the resulting model
pub fn get_atmosphere<P: AsRef<Path>>(path: P) -> Atmosphere {
    let mut file = File::open(path).unwrap();
    let mut contents = String::new();
    file.read_to_string(&mut contents).unwrap();

    let atm_def = parse_atmosphere(&contents).unwrap().1;
    Atmosphere::from_def(atm_def)
}

#[test]
fn test_find_t0() {
    let lapses = [0.0, 0.0065];
    let alts = [11e3, 20e3];
    let first_lapse = -0.0065;

    assert_eq!(
        Atmosphere::find_t0(0.0, 288.0, first_lapse, &lapses, &alts),
        288.0 - 11e3 * 0.0065
    );
    assert_eq!(
        Atmosphere::find_t0(12e3, 288.0, first_lapse, &lapses, &alts),
        288.0
    );
    assert_eq!(
        Atmosphere::find_t0(20e3, 288.0, first_lapse, &lapses, &alts),
        288.0
    );
    assert_eq!(
        Atmosphere::find_t0(28e3, 288.0, first_lapse, &lapses, &alts),
        288.0 - 0.0065 * 8e3
    );
}

#[test]
fn test_find_p0() {
    let lapses = [0.0, 0.0065];
    let alts = [11e3, 20e3];
    let first_lapse = -0.0065;
    let temps = [216.5, 216.5];

    assert_eq!(
        Atmosphere::find_p0(0.0, 1000.0, first_lapse, &lapses, &alts, &temps),
        223.15930527353143
    );
    assert_eq!(
        Atmosphere::find_p0(12e3, 1000.0, first_lapse, &lapses, &alts, &temps),
        1170.929298569465
    );
    assert_eq!(
        Atmosphere::find_p0(20e3, 1000.0, first_lapse, &lapses, &alts, &temps),
        4137.862509340013
    );
    assert_eq!(
        Atmosphere::find_p0(28e3, 1000.0, first_lapse, &lapses, &alts, &temps),
        12827.117117694357
    );
}

#[test]
fn test_us76() {
    let atmosphere = get_atmosphere("example_configs/us76.config");
    assert_eq!(atmosphere.pressure(0.0), 101325.0);
    assert_eq!(atmosphere.temperature(0.0), 288.0);
}

const US76: &str = "pressure(0) = 101325\
                    \
                    temperature:\
                    at(0) = 288\
                    lapse() = -0.0065\
                    lapse(11e3) = 0\
                    lapse(20e3) = 0.001\
                    lapse(32e3) = 0.0028\
                    lapse(47e3) = 0\
                    lapse(51e3) = -0.0028\
                    lapse(71e3) = -0.002\
                    lapse(84.852e3) = 0";

/// Returns the US-1976 standard model of the Earth's atmosphere.
///
/// The temperatures are expressed in kelvins (K), and the pressure in hectopascals (hPa).
pub fn us76_atmosphere() -> Atmosphere {
    let atm_def = parse_atmosphere(US76).unwrap().1;
    Atmosphere::from_def(atm_def)
}
