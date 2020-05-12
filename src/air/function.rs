use cubic_splines::Spline;

/// Something that can be evaluated for various values of an argument
pub trait Function {
    fn eval(&self, x: f64) -> f64;
    fn derivative(&self, x: f64) -> f64;
}

pub struct FunctionArray {
    interval_edges: Vec<f64>,
    // there should be one more function than interval edge;
    // the first function is on the interval (-infty, interval_edges[0]);
    // the last is on the interval (interval_edges.last(), infty)
    functions: Vec<Box<dyn Function>>,
}

impl Function for FunctionArray {
    fn eval(&self, x: f64) -> f64 {
        let index = self
            .interval_edges
            .binary_search_by(|probe| probe.partial_cmp(&x).unwrap());
        match index {
            Ok(i) => self.functions[i].eval(x),
            Err(i) => self.functions[i].eval(x),
        }
    }

    fn derivative(&self, x: f64) -> f64 {
        let index = self
            .interval_edges
            .binary_search_by(|probe| probe.partial_cmp(&x).unwrap());
        match index {
            Ok(i) => self.functions[i].derivative(x),
            Err(i) => self.functions[i].derivative(x),
        }
    }
}

impl FunctionArray {
    /// Creates a new FunctionArray
    /// Note: `edges` should yield f64s in ascending order
    pub fn new<I, F>(edges: I, functions: F) -> Self
    where
        I: IntoIterator<Item = f64>,
        F: IntoIterator<Item = Box<dyn Function>>,
    {
        let interval_edges: Vec<_> = edges.into_iter().collect();
        let functions: Vec<_> = functions.into_iter().collect();
        assert_eq!(functions.len(), interval_edges.len() + 1);
        Self {
            interval_edges,
            functions,
        }
    }
}

impl Function for Spline {
    fn eval(&self, x: f64) -> f64 {
        Spline::eval(self, x)
    }

    fn derivative(&self, x: f64) -> f64 {
        self.eval_derivative(x)
    }
}

pub struct Identity;

impl Function for Identity {
    fn eval(&self, x: f64) -> f64 {
        x
    }

    fn derivative(&self, _x: f64) -> f64 {
        1.0
    }
}

/// a * x + b
pub struct Linear {
    a: f64,
    b: f64,
    arg: Box<dyn Function>,
}

impl Linear {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            a,
            b,
            arg: Box::new(Identity),
        }
    }

    pub fn wrapping(a: f64, b: f64, arg: Box<dyn Function>) -> Self {
        Self { a, b, arg }
    }
}

impl Function for Linear {
    fn eval(&self, x: f64) -> f64 {
        self.a * self.arg.eval(x) + self.b
    }

    fn derivative(&self, x: f64) -> f64 {
        self.a * self.arg.derivative(x)
    }
}

/// a * exp(b * x)
pub struct Exponential {
    a: f64,
    b: f64,
    arg: Box<dyn Function>,
}

impl Exponential {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            a,
            b,
            arg: Box::new(Identity),
        }
    }

    pub fn wrapping(a: f64, b: f64, arg: Box<dyn Function>) -> Self {
        Self { a, b, arg }
    }
}

impl Function for Exponential {
    fn eval(&self, x: f64) -> f64 {
        self.a * (self.b * self.arg.eval(x)).exp()
    }

    fn derivative(&self, x: f64) -> f64 {
        self.eval(x) * self.b * self.arg.derivative(x)
    }
}

/// x^a
pub struct Power {
    a: f64,
    arg: Box<dyn Function>,
}

impl Power {
    pub fn new(a: f64) -> Self {
        Self {
            a,
            arg: Box::new(Identity),
        }
    }

    pub fn wrapping(a: f64, arg: Box<dyn Function>) -> Self {
        Self { a, arg }
    }
}

impl Function for Power {
    fn eval(&self, x: f64) -> f64 {
        self.arg.eval(x).powf(self.a)
    }

    fn derivative(&self, x: f64) -> f64 {
        self.a * self.arg.eval(x).powf(self.a - 1.0) * self.arg.derivative(x)
    }
}
