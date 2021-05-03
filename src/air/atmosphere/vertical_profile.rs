use cubic_splines::{BoundaryCondition, CubicPoly, Spline};

#[derive(Clone, Debug, PartialEq)]
enum VerticalFunction {
    /// T(h) = a*h + b
    Linear {
        a: f64,
        b: f64,
    },
    Cubic(CubicPoly<f64>),
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct VerticalProfile {
    altitude_interval_ends: Vec<f64>,
    interval_functions: Vec<VerticalFunction>,
}

#[derive(Clone, Debug)]
pub enum FunctionDef {
    Linear {
        gradient: f64,
    },
    Spline {
        points: Vec<(f64, f64)>,
        boundary_condition: BoundaryCondition<f64>,
    },
}

impl FunctionDef {
    fn into_intermediate(
        self,
        start_alt: Option<f64>,
        end_alt: Option<f64>,
    ) -> (Vec<f64>, Vec<IntermediateFunctionDef>) {
        match self {
            FunctionDef::Linear { gradient } => (
                start_alt.into_iter().collect(),
                vec![IntermediateFunctionDef::Linear {
                    gradient,
                    fixed_point: None,
                }],
            ),
            FunctionDef::Spline {
                points,
                boundary_condition,
            } => {
                let spline = Spline::new(points, boundary_condition);
                let mut alts = vec![];
                let mut funs = vec![];
                if start_alt.map_or(true, |start_alt| start_alt < spline.min_x()) {
                    if let Some(start_alt) = start_alt {
                        alts.push(start_alt);
                    }
                    funs.push(IntermediateFunctionDef::Linear {
                        gradient: spline.derivative_start(),
                        fixed_point: Some((spline.min_x(), spline.eval(spline.min_x()))),
                    });
                }
                for (mut start, end, poly) in spline.polynomials() {
                    if let Some(start_alt) = start_alt {
                        if start_alt > end {
                            continue;
                        }
                        if start_alt > start {
                            start = start_alt;
                        }
                    }
                    if let Some(end_alt) = end_alt {
                        if end_alt < start {
                            continue;
                        }
                    }
                    alts.push(start);
                    funs.push(IntermediateFunctionDef::Cubic { poly });
                }
                if end_alt.map_or(true, |end_alt| end_alt > spline.max_x()) {
                    alts.push(spline.max_x());
                    funs.push(IntermediateFunctionDef::Linear {
                        gradient: spline.derivative_end(),
                        fixed_point: Some((spline.max_x(), spline.eval(spline.max_x()))),
                    });
                }
                (alts, funs)
            }
        }
    }
}

#[derive(Clone, Debug)]
pub enum IntermediateFunctionDef {
    Linear {
        gradient: f64,
        fixed_point: Option<(f64, f64)>,
    },
    Cubic {
        poly: CubicPoly<f64>,
    },
}

impl IntermediateFunctionDef {
    fn get(&self, x: f64) -> Option<f64> {
        match self {
            IntermediateFunctionDef::Linear {
                gradient,
                fixed_point: Some((x0, y0)),
            } => Some(y0 + (x - x0) * gradient),
            IntermediateFunctionDef::Linear { .. } => None,
            IntermediateFunctionDef::Cubic { poly } => Some(poly.eval(x)),
        }
    }

    fn has_fixed_point(&self) -> bool {
        match self {
            IntermediateFunctionDef::Linear { fixed_point, .. } => fixed_point.is_some(),
            IntermediateFunctionDef::Cubic { .. } => true,
        }
    }

    fn into_function(self) -> VerticalFunction {
        match self {
            IntermediateFunctionDef::Linear {
                gradient,
                fixed_point,
            } => {
                let (x, y) = fixed_point.expect("should have a fixed point now");
                VerticalFunction::Linear {
                    a: gradient,
                    b: y - gradient * x,
                }
            }
            IntermediateFunctionDef::Cubic { poly } => VerticalFunction::Cubic(poly),
        }
    }
}

#[derive(Clone, Debug)]
pub struct VerticalProfileBuilder {
    interval_ends: Vec<f64>,
    function_defs: Vec<FunctionDef>,
    // value at a specific altitude - required if all the functions in the
    // definition are gradients
    fixed_value: Option<(f64, f64)>,
}

#[derive(Clone, Copy, Debug)]
pub enum FixedPoint {
    None,
    Spline(usize),
    FixedTemp,
}

impl VerticalProfileBuilder {
    pub fn new(first_fun: FunctionDef) -> Self {
        Self {
            interval_ends: vec![],
            function_defs: vec![first_fun],
            fixed_value: None,
        }
    }

    pub fn with_next_function(mut self, altitude: f64, function: FunctionDef) -> Self {
        self.interval_ends.push(altitude);
        self.function_defs.push(function);
        self
    }

    pub fn with_fixed_value(mut self, altitude: f64, val: f64) -> Self {
        self.fixed_value = Some((altitude, val));
        self
    }

    pub fn build(self) -> Result<VerticalProfile, VerticalProfileError> {
        let Self {
            interval_ends,
            function_defs,
            fixed_value,
        } = self;
        let (interval_ends, mut intermediate_function_defs) =
            Self::generate_intermediate_function_defs(interval_ends, function_defs);
        Self::fill_fixed_points(&interval_ends, &mut intermediate_function_defs, fixed_value)?;
        Ok(VerticalProfile {
            altitude_interval_ends: interval_ends,
            interval_functions: intermediate_function_defs
                .into_iter()
                .map(|fun_def| fun_def.into_function())
                .collect(),
        })
    }

    fn generate_intermediate_function_defs(
        mut interval_ends: Vec<f64>,
        mut function_defs: Vec<FunctionDef>,
    ) -> (Vec<f64>, Vec<IntermediateFunctionDef>) {
        let mut new_interval_ends = vec![];
        let mut new_function_defs = vec![];
        let first_function_def = function_defs.remove(0);
        let (alts, funs) =
            first_function_def.into_intermediate(None, interval_ends.get(0).cloned());
        new_interval_ends.extend(alts);
        new_function_defs.extend(funs);
        while !interval_ends.is_empty() {
            let start = interval_ends.remove(0);
            let fun = function_defs.remove(0);
            let (alts, funs) = fun.into_intermediate(Some(start), interval_ends.get(0).cloned());
            new_interval_ends.extend(alts);
            new_function_defs.extend(funs);
        }
        (new_interval_ends, new_function_defs)
    }

    fn fill_fixed_points(
        interval_ends: &[f64],
        function_defs: &mut [IntermediateFunctionDef],
        fixed_value: Option<(f64, f64)>,
    ) -> Result<(), VerticalProfileError> {
        for index in 0..function_defs.len() {
            Self::fill_fixed_point(interval_ends, function_defs, index, fixed_value)?;
        }
        Ok(())
    }

    fn fill_fixed_point(
        interval_ends: &[f64],
        function_defs: &mut [IntermediateFunctionDef],
        index: usize,
        fixed_value: Option<(f64, f64)>,
    ) -> Result<(), VerticalProfileError> {
        const EPSILON: f64 = 1e-4;

        let has_fixed_point_below = index.checked_sub(1).map_or(false, |index_below| {
            function_defs[index_below].has_fixed_point()
        });
        let has_fixed_point_above =
            (index + 1 < function_defs.len()) && function_defs[index + 1].has_fixed_point();

        // if there was a fixed_value specified and it's in our interval, set it as our
        // fixed point
        if let (Some((x, _)), IntermediateFunctionDef::Linear { fixed_point, .. }) =
            (fixed_value, &mut function_defs[index])
        {
            if index
                .checked_sub(1)
                .map_or(true, |ib| interval_ends[ib] <= x)
                && (index >= interval_ends.len() || interval_ends[index] >= x)
            {
                *fixed_point = fixed_value;
            }
        }

        let has_fixed_point = function_defs[index].has_fixed_point();
        if !has_fixed_point_below && !has_fixed_point_above && !has_fixed_point {
            if index + 1 < function_defs.len() {
                Self::fill_fixed_point(interval_ends, function_defs, index + 1, fixed_value)?;
            } else {
                return Err(VerticalProfileError::NoFixedPoint);
            }
        }

        let point_below = index.checked_sub(1).and_then(|index_below| {
            let x = interval_ends[index_below];
            function_defs[index_below].get(x).map(|y| (x, y))
        });
        let point_above = if index + 1 < function_defs.len() {
            let x = interval_ends[index];
            function_defs[index + 1].get(x).map(|y| (x, y))
        } else {
            None
        };

        match &mut function_defs[index] {
            IntermediateFunctionDef::Linear {
                gradient,
                fixed_point,
            } => {
                // check the consistency of the previous function's fixed point and ours
                if let (Some((x1, y1)), Some((x2, y2))) = (point_below, *fixed_point) {
                    if ((y2 - y1) - *gradient * (x2 - x1)).abs() > EPSILON {
                        return Err(VerticalProfileError::FixedPointConflict {
                            index1: index - 1,
                            index2: index,
                            point1: (x1, y1),
                            point2: (x2, y2),
                            gradient: Some(*gradient),
                        });
                    }
                }
                // check the consistency of the next function's fixed point and ours
                if let (Some((x1, y1)), Some((x2, y2))) = (*fixed_point, point_above) {
                    if ((y2 - y1) - *gradient * (x2 - x1)).abs() > EPSILON {
                        return Err(VerticalProfileError::FixedPointConflict {
                            index1: index,
                            index2: index + 1,
                            point1: (x1, y1),
                            point2: (x2, y2),
                            gradient: Some(*gradient),
                        });
                    }
                }
                if fixed_point.is_none() {
                    if let Some(point) = point_below {
                        *fixed_point = Some(point);
                        Ok(())
                    } else if let Some(point) = point_above {
                        *fixed_point = Some(point);
                        Ok(())
                    } else {
                        Err(VerticalProfileError::NoFixedPoint)
                    }
                } else {
                    Ok(())
                }
            }
            IntermediateFunctionDef::Cubic { poly } => {
                // if the fixed value is in our interval, check its consistency with the polynomial
                if let Some((x, y)) = fixed_value {
                    if index
                        .checked_sub(1)
                        .map_or(true, |ib| interval_ends[ib] <= x)
                        && (index >= interval_ends.len() || interval_ends[index] >= x)
                        && (y - poly.eval(x)).abs() > EPSILON
                    {
                        return Err(VerticalProfileError::FixedPointConflict {
                            index1: index,
                            index2: index,
                            point1: (x, y),
                            point2: (x, poly.eval(x)),
                            gradient: None,
                        });
                    }
                }
                // check the consistency of the previous function's fixed point and ours
                if let Some((x, y)) = point_below {
                    if (y - poly.eval(x)).abs() > EPSILON {
                        return Err(VerticalProfileError::FixedPointConflict {
                            index1: index - 1,
                            index2: index,
                            point1: (x, y),
                            point2: (x, poly.eval(x)),
                            gradient: None,
                        });
                    }
                }
                // check the consistency of the next function's fixed point and ours
                if let Some((x, y)) = point_above {
                    if (y - poly.eval(x)).abs() > EPSILON {
                        return Err(VerticalProfileError::FixedPointConflict {
                            index1: index,
                            index2: index + 1,
                            point1: (x, poly.eval(x)),
                            point2: (x, y),
                            gradient: None,
                        });
                    }
                }
                Ok(())
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VerticalProfileError {
    NoFixedPoint,
    FixedPointConflict {
        index1: usize,
        index2: usize,
        point1: (f64, f64),
        point2: (f64, f64),
        gradient: Option<f64>,
    },
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn should_build_correctly_with_fixed_point_in_first_interval() {
        let _ = VerticalProfileBuilder::new(FunctionDef::Linear { gradient: -0.0065 })
            .with_next_function(11e3, FunctionDef::Linear { gradient: 0.0 })
            .with_next_function(15e3, FunctionDef::Linear { gradient: 0.0065 })
            .with_fixed_value(0.0, 288.0)
            .build()
            .expect("should build correctly");
    }

    #[test]
    fn should_build_correctly_with_fixed_point_in_last_interval() {
        let _ = VerticalProfileBuilder::new(FunctionDef::Linear { gradient: -0.0065 })
            .with_next_function(11e3, FunctionDef::Linear { gradient: 0.0 })
            .with_next_function(15e3, FunctionDef::Linear { gradient: 0.0065 })
            .with_fixed_value(16e3, 218.0)
            .build()
            .expect("should build correctly");
    }

    #[test]
    fn should_build_correctly_with_fixed_point_in_the_middle() {
        let _ = VerticalProfileBuilder::new(FunctionDef::Linear { gradient: -0.0065 })
            .with_next_function(11e3, FunctionDef::Linear { gradient: 0.0 })
            .with_next_function(15e3, FunctionDef::Linear { gradient: 0.0065 })
            .with_fixed_value(12e3, 218.0)
            .build()
            .expect("should build correctly");
    }

    #[test]
    fn should_build_correctly_with_only_spline() {
        let _ = VerticalProfileBuilder::new(FunctionDef::Spline {
            points: vec![(0.0, 0.0), (10.0, -2.0), (15.0, 3.0)],
            boundary_condition: BoundaryCondition::Natural,
        })
        .build()
        .expect("should build correctly");
    }

    #[test]
    fn should_build_correctly_with_spline_with_fixed_point() {
        let _ = VerticalProfileBuilder::new(FunctionDef::Spline {
            points: vec![(0.0, 0.0), (10.0, -2.0), (15.0, 3.0)],
            boundary_condition: BoundaryCondition::Natural,
        })
        .with_fixed_value(10.0, -2.0)
        .build()
        .expect("should build correctly");
    }

    #[test]
    fn should_build_correctly_with_spline_and_linear() {
        let _ = VerticalProfileBuilder::new(FunctionDef::Spline {
            points: vec![(0.0, 0.0), (10.0, -2.0), (15.0, 3.0)],
            boundary_condition: BoundaryCondition::Natural,
        })
        .with_next_function(16.0, FunctionDef::Linear { gradient: 3.0 })
        .build()
        .expect("should build correctly");
    }

    #[test]
    fn should_build_correctly_with_linear_and_spline() {
        let _ = VerticalProfileBuilder::new(FunctionDef::Linear { gradient: 3.0 })
            .with_next_function(
                -1.0,
                FunctionDef::Spline {
                    points: vec![(0.0, 0.0), (10.0, -2.0), (15.0, 3.0)],
                    boundary_condition: BoundaryCondition::Natural,
                },
            )
            .build()
            .expect("should build correctly");
    }

    #[test]
    fn should_fail_if_linear_without_fixed_value() {
        let result = VerticalProfileBuilder::new(FunctionDef::Linear { gradient: 3.1 })
            .with_next_function(0.0, FunctionDef::Linear { gradient: 3.0 })
            .build();
        assert_eq!(result, Err(VerticalProfileError::NoFixedPoint));
    }

    #[test]
    fn should_fail_if_conflict() {
        let result = VerticalProfileBuilder::new(FunctionDef::Linear { gradient: 3.0 })
            .with_next_function(
                0.0,
                FunctionDef::Spline {
                    points: vec![(0.0, 0.0), (10.0, -2.0), (15.0, 3.0)],
                    boundary_condition: BoundaryCondition::Natural,
                },
            )
            .with_fixed_value(-2.0, 0.0)
            .build();
        assert_eq!(
            result,
            Err(VerticalProfileError::FixedPointConflict {
                index1: 0,
                index2: 1,
                point1: (-2.0, 0.0),
                point2: (0.0, 0.0),
                gradient: Some(3.0),
            })
        );
    }
}
