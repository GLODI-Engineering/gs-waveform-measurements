//! The generic `(time, value)` time series every measurement in this crate operates on, plus
//! the small set of primitives (linear interpolation, windowing) that all of them share.
//!
//! Every measurement in this crate treats its input as **piecewise-linear (PWL) between real
//! samples** — the same assumption a variable-timestep transient solver's own output naturally
//! satisfies, since the solver only stores the solution at the (non-uniformly spaced) instants
//! it actually solved for. No measurement here resamples to a uniform grid; every quantity is
//! computed by walking the real samples and, where a threshold or a window boundary falls
//! strictly between two samples, linearly interpolating between exactly those two.

use core::fmt;

/// A single `(time, value)` observation of a signal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sample {
    /// The independent variable — simulation time, in this crate's intended use, but nothing
    /// here assumes that; any monotonically increasing coordinate works.
    pub t: f64,
    /// The signal's value at `t`.
    pub v: f64,
}

impl Sample {
    /// Constructs a new sample.
    pub fn new(t: f64, v: f64) -> Self {
        Self { t, v }
    }
}

/// A measurement window, mirroring ngspice/Xyce's `FROM=`/`TO=`/`TD=` qualifiers.
///
/// The effective window is the intersection `[max(FROM, TD), TO]` of the series' own domain
/// with these bounds (Xyce's own documented intent: "the measurement window \[is\] the
/// intersection of the FROM-TO and TD windows, if both are specified"). Any bound left as
/// `None` is unconstrained on that side.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Window {
    /// Start of the window (`FROM=`). `None` means unconstrained (the series' own start).
    pub from: Option<f64>,
    /// End of the window (`TO=`). `None` means unconstrained (the series' own end).
    pub to: Option<f64>,
    /// A delay before which the measurement is not evaluated (`TD=`). Combines with `from` by
    /// taking the later of the two, per the Xyce semantics above.
    pub td: Option<f64>,
}

impl Window {
    /// The unconstrained window: the entire series.
    pub const FULL: Window = Window {
        from: None,
        to: None,
        td: None,
    };

    /// A window starting at `from` and ending at `to` (both bounds inclusive), with no `TD`.
    pub fn new(from: f64, to: f64) -> Self {
        Self {
            from: Some(from),
            to: Some(to),
            td: None,
        }
    }

    /// The effective lower bound: `max(from, td)`, treating an absent bound as unconstrained.
    fn effective_from(&self) -> f64 {
        match (self.from, self.td) {
            (Some(a), Some(b)) => a.max(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => f64::NEG_INFINITY,
        }
    }

    fn effective_to(&self) -> f64 {
        self.to.unwrap_or(f64::INFINITY)
    }
}

/// Every fallible operation in this crate reports one of these reasons rather than silently
/// substituting a default value (unlike ngspice/Xyce's `DEFAULT_VAL=` semantics, which return a
/// user-chosen fallback and log a failure) — callers that want Xyce's "return a default on
/// failure" behavior can layer `.unwrap_or(default)` on top of a `Result`-returning call here.
#[derive(Debug, Clone, PartialEq)]
pub enum MeasureError {
    /// The input series had fewer than two samples, so no interval-based quantity is defined.
    EmptySeries,
    /// The requested window (after intersecting `FROM`/`TO`/`TD` with the series' own domain)
    /// contains fewer than two samples, so no interval-based quantity is defined over it.
    EmptyWindow,
    /// A crossing, threshold, or fixed time was requested but never occurs in the (windowed)
    /// series.
    NotFound,
    /// A parameter combination is invalid independent of the data (e.g. zero harmonics
    /// requested, or a non-monotonic input series).
    InvalidParameter(&'static str),
}

impl fmt::Display for MeasureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MeasureError::EmptySeries => write!(f, "series has fewer than two samples"),
            MeasureError::EmptyWindow => {
                write!(f, "measurement window contains fewer than two samples")
            }
            MeasureError::NotFound => write!(f, "no matching crossing/event found in the window"),
            MeasureError::InvalidParameter(msg) => write!(f, "invalid parameter: {msg}"),
        }
    }
}

impl std::error::Error for MeasureError {}

/// Convenience alias used throughout this crate.
pub type Result<T> = core::result::Result<T, MeasureError>;

/// Linearly interpolates `samples` (which must be sorted by non-decreasing `t`) at time `t`.
///
/// Returns `None` if `t` falls outside `[samples[0].t, samples[last].t]`. This is the single
/// primitive every "value at an arbitrary time" measurement (`FIND ... AT=`, window-boundary
/// clipping, crossing interpolation) is built from.
pub fn interpolate(samples: &[Sample], t: f64) -> Option<f64> {
    if samples.is_empty() {
        return None;
    }
    let first = samples[0];
    let last = samples[samples.len() - 1];
    if t < first.t || t > last.t {
        return None;
    }
    if t == first.t {
        return Some(first.v);
    }
    if t == last.t {
        return Some(last.v);
    }
    // Binary search for the segment straddling t.
    let idx = match samples.binary_search_by(|s| s.t.partial_cmp(&t).unwrap()) {
        Ok(i) => return Some(samples[i].v),
        Err(i) => i,
    };
    let a = samples[idx - 1];
    let b = samples[idx];
    let frac = (t - a.t) / (b.t - a.t);
    Some(a.v + frac * (b.v - a.v))
}

/// Verifies `samples` is sorted by strictly increasing `t`, as every measurement in this crate
/// requires.
pub fn check_monotonic(samples: &[Sample]) -> Result<()> {
    if samples.windows(2).any(|w| w[1].t <= w[0].t) {
        return Err(MeasureError::InvalidParameter(
            "samples must be sorted by strictly increasing t",
        ));
    }
    Ok(())
}

/// Clips `samples` to the effective window, inserting an interpolated `Sample` at each boundary
/// that falls strictly between two real samples, so that trapezoidal/analytic integration over
/// the returned slice is exact up to the window edges rather than snapping to the nearest real
/// sample.
///
/// Returns [`MeasureError::EmptyWindow`] if fewer than two samples remain (including the
/// interpolated boundary points).
pub fn clip(samples: &[Sample], window: &Window) -> Result<Vec<Sample>> {
    if samples.is_empty() {
        return Err(MeasureError::EmptySeries);
    }
    check_monotonic(samples)?;

    let lo = window.effective_from().max(samples[0].t);
    let hi = window.effective_to().min(samples[samples.len() - 1].t);
    if lo >= hi {
        return Err(MeasureError::EmptyWindow);
    }

    let mut out = Vec::with_capacity(samples.len());
    if let Some(v) = interpolate(samples, lo) {
        out.push(Sample::new(lo, v));
    }
    for s in samples {
        if s.t > lo && s.t < hi {
            out.push(*s);
        }
    }
    if let Some(v) = interpolate(samples, hi) {
        out.push(Sample::new(hi, v));
    }

    if out.len() < 2 {
        return Err(MeasureError::EmptyWindow);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolate_midpoint() {
        let s = [Sample::new(0.0, 0.0), Sample::new(2.0, 10.0)];
        assert_eq!(interpolate(&s, 1.0), Some(5.0));
    }

    #[test]
    fn interpolate_outside_domain_is_none() {
        let s = [Sample::new(0.0, 0.0), Sample::new(2.0, 10.0)];
        assert_eq!(interpolate(&s, -1.0), None);
        assert_eq!(interpolate(&s, 3.0), None);
    }

    #[test]
    fn clip_inserts_interpolated_boundaries() {
        let s = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 10.0),
            Sample::new(2.0, 0.0),
        ];
        let w = Window::new(0.5, 1.5);
        let c = clip(&s, &w).unwrap();
        assert_eq!(c[0], Sample::new(0.5, 5.0));
        assert_eq!(c[1], Sample::new(1.0, 10.0));
        assert_eq!(c[2], Sample::new(1.5, 5.0));
    }
}
