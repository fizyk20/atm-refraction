use nom::{self, types::CompleteStr};

named!(float <CompleteStr, f64>, map!(
        re_find!(r"^(\+|-)?\d+(\.\d+)?(e(\+|-)?\d+)?"),
        |s| s.parse().unwrap()));

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PressureDef {
    pub start_h: f64,
    pub start_p: f64,
}

named!(pressure_def <CompleteStr, PressureDef>, ws!(do_parse!(
        tag!("pressure") >>
        char!('(') >>
        start_h: float >>
        char!(')') >>
        tag!("=") >>
        start_p: float >>
        (PressureDef { start_h, start_p })
        )));

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LapseDef {
    pub start_h: Option<f64>,
    pub lapse: f64,
}

named!(lapse_def <CompleteStr, LapseDef>, ws!(do_parse!(
            tag!("lapse") >>
            char!('(') >>
            start_h: opt!(float) >>
            char!(')') >>
            tag!("=") >>
            lapse: float >>
            (LapseDef { start_h, lapse })
    )));

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TemperatureAt {
    pub start_h: f64,
    pub start_t: f64,
}

named!(temperature_at <CompleteStr, TemperatureAt>, ws!(do_parse!(
        tag!("at") >>
        char!('(') >>
        start_h: float >>
        char!(')') >>
        tag!("=") >>
        start_t: float >>
        (TemperatureAt { start_h, start_t })
        )));

#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureDef {
    pub start: TemperatureAt,
    pub lapses: Vec<LapseDef>,
}

named!(temperature_def <CompleteStr, TemperatureDef>, ws!(do_parse!(
            tag!("temperature:") >>
            start: temperature_at >>
            lapses: many1!(lapse_def) >>
            (TemperatureDef { start, lapses })
            )));

#[derive(Debug, Clone, PartialEq)]
pub struct AtmosphereDef {
    pub temperature: TemperatureDef,
    pub pressure: PressureDef,
}

named!(atmosphere <CompleteStr, AtmosphereDef>, ws!(do_parse!(
            pressure: pressure_def >>
            temperature: temperature_def >>
            (AtmosphereDef { temperature, pressure })
        )));

pub fn parse_atmosphere(txt: &str) -> nom::IResult<CompleteStr, AtmosphereDef> {
    atmosphere(CompleteStr(txt))
}

#[test]
fn test_float() {
    assert_eq!(float(CompleteStr("1")), Ok((CompleteStr(""), 1.0)));
    assert_eq!(float(CompleteStr("-3")), Ok((CompleteStr(""), -3.0)));
    assert_eq!(float(CompleteStr("+50")), Ok((CompleteStr(""), 50.0)));
    assert_eq!(float(CompleteStr("0.3")), Ok((CompleteStr(""), 0.3)));
    assert_eq!(float(CompleteStr("1.8e3")), Ok((CompleteStr(""), 1.8e3)));
    assert_eq!(float(CompleteStr("-5.03")), Ok((CompleteStr(""), -5.03)));
    assert_eq!(
        float(CompleteStr("-5.03e-2")),
        Ok((CompleteStr(""), -5.03e-2))
    );
}

#[test]
fn test_pressure_def() {
    assert_eq!(
        pressure_def(CompleteStr("pressure(0) = 101325.3")),
        Ok((
            CompleteStr(""),
            PressureDef {
                start_h: 0.0,
                start_p: 101325.3
            }
        ))
    );
    assert_eq!(
        pressure_def(CompleteStr("pressure (0) = 101325")),
        Ok((
            CompleteStr(""),
            PressureDef {
                start_h: 0.0,
                start_p: 101325.0
            }
        ))
    );
    assert_eq!(
        pressure_def(CompleteStr("pressure( 0 ) = 101325")),
        Ok((
            CompleteStr(""),
            PressureDef {
                start_h: 0.0,
                start_p: 101325.0
            }
        ))
    );
    assert_eq!(
        pressure_def(CompleteStr("pressure(0)=101325")),
        Ok((
            CompleteStr(""),
            PressureDef {
                start_h: 0.0,
                start_p: 101325.0
            }
        ))
    );
}

#[test]
fn test_lapse_def() {
    assert_eq!(
        lapse_def(CompleteStr("lapse() = 0.32")),
        Ok((
            CompleteStr(""),
            LapseDef {
                start_h: None,
                lapse: 0.32
            }
        ))
    );
    assert_eq!(
        lapse_def(CompleteStr("     lapse(11e3) = 0.0")),
        Ok((
            CompleteStr(""),
            LapseDef {
                start_h: Some(11e3),
                lapse: 0.0
            }
        ))
    );
}

#[test]
fn test_temperature_def() {
    assert_eq!(
        temperature_def(CompleteStr(
            "temperature:\nat(0) = 288\nlapse() = -0.0065\nlapse(11e3) = 0"
        )),
        Ok((
            CompleteStr(""),
            TemperatureDef {
                start: TemperatureAt {
                    start_h: 0.0,
                    start_t: 288.0,
                },
                lapses: vec![
                    LapseDef {
                        start_h: None,
                        lapse: -0.0065
                    },
                    LapseDef {
                        start_h: Some(11e3),
                        lapse: 0.0
                    }
                ]
            }
        ))
    );
}
