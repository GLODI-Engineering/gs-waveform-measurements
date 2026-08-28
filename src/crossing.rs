//! Threshold-crossing measurements: `FIND ... AT=`, `WHEN`/`FIND ... WHEN`, `DERIV`, and
//! `TRIG`/`TARG`.
//!
//! Every crossing here is found by **linear interpolation strictly between the two adjacent
//! real samples that straddle it** — never by resampling the trace to a uniform grid. This
//! mirrors ngspice/Xyce's own documented "level-crossing" mode: walking consecutive samples and
//! testing `(current - crossVal)` against `(previous - crossVal)` for a sign change.

use crate::sample::{check_monotonic, clip, interpolate, MeasureError, Result, Sample, Window};
use crate::scalar;

/// Which direction(s) of threshold crossing to accept — ngspice/Xyce's `RISE=`, `FALL=`,
/// `CROSS=` qualifiers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Edge {
    /// Only rising crossings (signal goes from below the threshold to at/above it) — `RISE=`.
    Rising,
    /// Only falling crossings (signal goes from above the threshold to at/below it) — `FALL=`.
    Falling,
    /// Either direction — `CROSS=`.
    Either,
}

/// Which occurrence of a matching crossing to use — ngspice/Xyce's `RISE=n`/`FALL=n`/`CROSS=n`
/// (1-based count) versus `RISE=LAST`/`FALL=LAST`/`CROSS=LAST`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Occurrence {
    /// The `n`-th matching crossing, 1-based (`RISE=1` is the *first* rising crossing).
    Nth(usize),
    /// The last matching crossing found in the window.
    Last,
}

/// The right-hand side of a `WHEN <variable> = <value>|<variable2>` clause.
pub enum Threshold<'a> {
    /// A fixed threshold value.
    Value(f64),
    /// Another time series; the crossing is where the two series' difference changes sign.
    Variable(&'a [Sample]),
}

/// `FIND <variable> AT=<time>`: the value of `samples` at an arbitrary time, by linear
/// interpolation between the two real samples straddling it.
pub fn find_at(samples: &[Sample], t: f64) -> Result<f64> {
    check_monotonic(samples)?;
    interpolate(samples, t).ok_or(MeasureError::NotFound)
}

/// `DERIV <variable> AT=<time>`: the derivative of `samples` at `t`.
///
/// Because the data is treated as piecewise-linear, the derivative is *constant* within the
/// segment containing `t` — it is exactly that segment's slope, not merely a finite-difference
/// approximation. If `t` lands exactly on an interior sample (a segment boundary), the natural
/// generalization is the central difference: the average of the slopes of the segment ending at
/// `t` and the segment starting at `t`.
pub fn deriv_at(samples: &[Sample], t: f64) -> Result<f64> {
    check_monotonic(samples)?;
    if samples.len() < 2 {
        return Err(MeasureError::EmptySeries);
    }
    if t < samples[0].t || t > samples[samples.len() - 1].t {
        return Err(MeasureError::NotFound);
    }
    let slope = |a: Sample, b: Sample| (b.v - a.v) / (b.t - a.t);

    // Find i such that samples[i].t <= t <= samples[i+1].t.
    let idx = match samples.binary_search_by(|s| s.t.partial_cmp(&t).unwrap()) {
        Ok(i) => i,
        Err(i) => i - 1,
    };
    let on_boundary = samples[idx].t == t;
    if on_boundary {
        let left = if idx > 0 {
            Some(slope(samples[idx - 1], samples[idx]))
        } else {
            None
        };
        let right = if idx + 1 < samples.len() {
            Some(slope(samples[idx], samples[idx + 1]))
        } else {
            None
        };
        return match (left, right) {
            (Some(l), Some(r)) => Ok(0.5 * (l + r)),
            (Some(l), None) => Ok(l),
            (None, Some(r)) => Ok(r),
            (None, None) => Err(MeasureError::EmptySeries),
        };
    }
    Ok(slope(samples[idx], samples[idx + 1]))
}

/// A single located crossing: the (interpolated) time it occurs at, and its direction.
#[derive(Debug, Clone, Copy)]
struct CrossingEvent {
    t: f64,
    edge: Edge,
}

/// Builds the "difference from threshold" series on the union of both inputs' sample times
/// (restricted to their overlapping domain), so a crossing between two series with different
/// sample grids is still found by interpolating each onto the other's real sample times, not by
/// resampling either to a synthetic uniform grid.
fn diff_series(var: &[Sample], threshold: &Threshold<'_>) -> Result<Vec<Sample>> {
    check_monotonic(var)?;
    match threshold {
        Threshold::Value(v) => {
            if var.is_empty() {
                return Err(MeasureError::EmptySeries);
            }
            Ok(var.iter().map(|s| Sample::new(s.t, s.v - v)).collect())
        }
        Threshold::Variable(other) => {
            check_monotonic(other)?;
            if var.is_empty() || other.is_empty() {
                return Err(MeasureError::EmptySeries);
            }
            let lo = var[0].t.max(other[0].t);
            let hi = var[var.len() - 1].t.min(other[other.len() - 1].t);
            if lo >= hi {
                return Err(MeasureError::EmptyWindow);
            }
            let mut times: Vec<f64> = var
                .iter()
                .chain(other.iter())
                .map(|s| s.t)
                .filter(|&t| t >= lo && t <= hi)
                .collect();
            times.sort_by(|a, b| a.partial_cmp(b).unwrap());
            times.dedup_by(|a, b| (*a - *b).abs() < 1e-15);
            let diffs: Vec<Sample> = times
                .into_iter()
                .filter_map(|t| {
                    let a = interpolate(var, t)?;
                    let b = interpolate(other, t)?;
                    Some(Sample::new(t, a - b))
                })
                .collect();
            if diffs.len() < 2 {
                return Err(MeasureError::EmptyWindow);
            }
            Ok(diffs)
        }
    }
}

/// Walks a difference series and reports every zero-crossing, using the same level-crossing
/// pseudocode ngspice/Xyce document: rising if `current >= 0 && previous < 0`, falling if
/// `current <= 0 && previous > 0`.
fn find_crossings(diff: &[Sample]) -> Vec<CrossingEvent> {
    let mut out = Vec::new();
    for w in diff.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (d0, d1) = (a.v, b.v);
        if d0 == 0.0 && d1 == 0.0 {
            continue;
        }
        let interp_t = |d0: f64, d1: f64| -> f64 {
            if d1 == d0 {
                a.t
            } else {
                a.t + (0.0 - d0) * (b.t - a.t) / (d1 - d0)
            }
        };
        if d1 >= 0.0 && d0 < 0.0 {
            out.push(CrossingEvent {
                t: interp_t(d0, d1),
                edge: Edge::Rising,
            });
        } else if d1 <= 0.0 && d0 > 0.0 {
            out.push(CrossingEvent {
                t: interp_t(d0, d1),
                edge: Edge::Falling,
            });
        }
    }
    out
}

fn select_crossing(events: &[CrossingEvent], edge: Edge, occurrence: Occurrence) -> Result<f64> {
    let matching: Vec<&CrossingEvent> = events
        .iter()
        .filter(|e| edge == Edge::Either || e.edge == edge)
        .collect();
    match occurrence {
        Occurrence::Last => matching.last().map(|e| e.t).ok_or(MeasureError::NotFound),
        Occurrence::Nth(n) => {
            if n == 0 {
                return Err(MeasureError::InvalidParameter(
                    "RISE/FALL/CROSS occurrence is 1-based; 0 is not valid",
                ));
            }
            matching
                .get(n - 1)
                .map(|e| e.t)
                .ok_or(MeasureError::NotFound)
        }
    }
}

/// `WHEN <variable> = <value>|<variable2>`: the *time* at which `var` crosses `threshold`,
/// selecting the requested occurrence and edge direction within `window`.
pub fn when(
    var: &[Sample],
    threshold: Threshold<'_>,
    window: &Window,
    edge: Edge,
    occurrence: Occurrence,
) -> Result<f64> {
    let diff = diff_series(var, &threshold)?;
    let diff = clip(&diff, window)?;
    let events = find_crossings(&diff);
    select_crossing(&events, edge, occurrence)
}

/// `FIND <find_var> WHEN <when_var> = <value>|<variable2>`: the value of `find_var` at the
/// moment `when_var` crosses `threshold`.
pub fn find_when(
    find_var: &[Sample],
    when_var: &[Sample],
    threshold: Threshold<'_>,
    window: &Window,
    edge: Edge,
    occurrence: Occurrence,
) -> Result<f64> {
    let t = when(when_var, threshold, window, edge, occurrence)?;
    interpolate(find_var, t).ok_or(MeasureError::NotFound)
}

/// `DERIV <variable> WHEN <when_var> = <value>|<variable2>`: the derivative of `variable` at the
/// moment `when_var` crosses `threshold`.
pub fn deriv_when(
    variable: &[Sample],
    when_var: &[Sample],
    threshold: Threshold<'_>,
    window: &Window,
    edge: Edge,
    occurrence: Occurrence,
) -> Result<f64> {
    let t = when(when_var, threshold, window, edge, occurrence)?;
    deriv_at(variable, t)
}

/// The trigger/target event specification shared by `TRIG` and `TARG` clauses.
pub enum EventSpec<'a> {
    /// `TRIG AT=<value>` / `TARG AT=<value>`: a fixed, already-known time (no search).
    At(f64),
    /// `TRIG <var1>=<var2>|<value>` / `TARG <var3>=<var4>|<value>`: a threshold crossing on
    /// `var`, with its own `RISE=`/`FALL=`/`CROSS=` selection and window.
    Crossing {
        /// The series being searched for the crossing.
        var: &'a [Sample],
        /// What `var` is being compared against.
        threshold: Threshold<'a>,
        /// The measurement window for this side's crossing search.
        window: Window,
        /// Which edge direction to accept.
        edge: Edge,
        /// Which occurrence to select.
        occurrence: Occurrence,
    },
    /// `TRIG <var> FRAC_MAX=<value>` / `TARG <var> FRAC_MAX=<value>`: a crossing of `frac *
    /// max(var)` (the peak computed over `window`) rather than an absolute level — useful when
    /// the signal's peak isn't known in advance (e.g. ensemble/rise-time measurements).
    FracMax {
        /// The series being searched.
        var: &'a [Sample],
        /// Fraction of `var`'s own peak value (over `window`) to use as the crossing level.
        frac: f64,
        /// The measurement window (both for computing the peak and for the crossing search).
        window: Window,
        /// Which edge direction to accept.
        edge: Edge,
        /// Which occurrence to select.
        occurrence: Occurrence,
    },
}

/// Resolves an [`EventSpec`] to a concrete time.
pub fn resolve_event(spec: &EventSpec<'_>) -> Result<f64> {
    match spec {
        EventSpec::At(t) => Ok(*t),
        EventSpec::Crossing {
            var,
            threshold,
            window,
            edge,
            occurrence,
        } => {
            let diff = diff_series(var, threshold)?;
            let diff = clip(&diff, window)?;
            let events = find_crossings(&diff);
            select_crossing(&events, *edge, *occurrence)
        }
        EventSpec::FracMax {
            var,
            frac,
            window,
            edge,
            occurrence,
        } => {
            let peak = scalar::max(var, window)?.value;
            let level = frac * peak;
            let diff = diff_series(var, &Threshold::Value(level))?;
            let diff = clip(&diff, window)?;
            let events = find_crossings(&diff);
            select_crossing(&events, *edge, *occurrence)
        }
    }
}

/// `TRIG ... TARG ...`: the time difference between a target event and a trigger event,
/// `t_targ - t_trig`.
pub fn trig_targ(trig: &EventSpec<'_>, targ: &EventSpec<'_>) -> Result<f64> {
    let t_trig = resolve_event(trig)?;
    let t_targ = resolve_event(targ)?;
    Ok(t_targ - t_trig)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A ramp from (0,0) to (10,10): WHEN v=4 must land exactly at t=4, by hand (v=t here).
    #[test]
    fn when_on_linear_ramp_is_exact() {
        let s = [Sample::new(0.0, 0.0), Sample::new(10.0, 10.0)];
        let t = when(
            &s,
            Threshold::Value(4.0),
            &Window::FULL,
            Edge::Rising,
            Occurrence::Nth(1),
        )
        .unwrap();
        assert!((t - 4.0).abs() < 1e-12, "got {t}");
    }

    /// Non-uniform Δt around the crossing: samples at 0, 0.1, 0.15, 1.0 with values 0, 1, 1.5,
    /// 10 (piecewise linear, slopes 10, 10, ~9.55). Crossing of 1.2 must fall in the segment
    /// [0.15,1.0]: hand-solved t = 0.15 + (1.2-1.5)/(10-1.5)*(1.0-0.15) = 0.15 - 0.3/8.5*0.85
    /// = 0.15 - 0.03 = 0.12.
    #[test]
    fn when_respects_nonuniform_spacing() {
        let s = [
            Sample::new(0.0, 0.0),
            Sample::new(0.1, 1.0),
            Sample::new(0.15, 1.5),
            Sample::new(1.0, 10.0),
        ];
        let t = when(
            &s,
            Threshold::Value(1.2),
            &Window::FULL,
            Edge::Rising,
            Occurrence::Nth(1),
        )
        .unwrap();
        let expected = 0.15 + (1.2 - 1.5) / (10.0 - 1.5) * (1.0 - 0.15);
        assert!((t - expected).abs() < 1e-12, "got {t}, expected {expected}");
    }

    /// A triangle wave 0 -> 1 (t=0..1) -> 0 (t=1..2): two rising... actually one rise, one fall.
    /// RISE=1 must hit the up-slope crossing of 0.5 at t=0.5; FALL=1 hits the down-slope
    /// crossing at t=1.5, by hand (symmetric triangle).
    #[test]
    fn rise_and_fall_selection() {
        let s = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 1.0),
            Sample::new(2.0, 0.0),
        ];
        let rise = when(
            &s,
            Threshold::Value(0.5),
            &Window::FULL,
            Edge::Rising,
            Occurrence::Nth(1),
        )
        .unwrap();
        assert!((rise - 0.5).abs() < 1e-12);
        let fall = when(
            &s,
            Threshold::Value(0.5),
            &Window::FULL,
            Edge::Falling,
            Occurrence::Nth(1),
        )
        .unwrap();
        assert!((fall - 1.5).abs() < 1e-12);
    }

    /// find_at on a ramp v=2t: at t=3, value must be 6, by hand.
    #[test]
    fn find_at_linear() {
        let s = [Sample::new(0.0, 0.0), Sample::new(10.0, 20.0)];
        assert!((find_at(&s, 3.0).unwrap() - 6.0).abs() < 1e-12);
    }

    /// deriv_at strictly inside a segment of slope 3 must return exactly 3 (PWL derivative is
    /// exact, not a finite-difference approximation).
    #[test]
    fn deriv_at_interior_is_exact_segment_slope() {
        let s = [Sample::new(0.0, 0.0), Sample::new(2.0, 6.0)];
        assert!((deriv_at(&s, 1.0).unwrap() - 3.0).abs() < 1e-12);
    }

    /// deriv_at on a boundary between a slope-1 segment and a slope-3 segment must return the
    /// central difference (1+3)/2 = 2, by hand.
    #[test]
    fn deriv_at_boundary_is_central_difference() {
        let s = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 1.0),
            Sample::new(2.0, 4.0),
        ];
        assert!((deriv_at(&s, 1.0).unwrap() - 2.0).abs() < 1e-12);
    }

    /// TRIG at t=1 (fixed), TARG at the crossing of a ramp v=t through 5 (t=5): result = 5-1=4.
    #[test]
    fn trig_targ_fixed_and_crossing() {
        let s = [Sample::new(0.0, 0.0), Sample::new(10.0, 10.0)];
        let trig = EventSpec::At(1.0);
        let targ = EventSpec::Crossing {
            var: &s,
            threshold: Threshold::Value(5.0),
            window: Window::FULL,
            edge: Edge::Rising,
            occurrence: Occurrence::Nth(1),
        };
        let dt = trig_targ(&trig, &targ).unwrap();
        assert!((dt - 4.0).abs() < 1e-12, "got {dt}");
    }

    /// FRAC_MAX rise time: a ramp 0->10 over [0,10]; 10%-90% rise time is (9-1)=8 by hand
    /// (crossings at t=1 and t=9 for a v=t ramp with peak 10).
    #[test]
    fn frac_max_rise_time() {
        let s = [Sample::new(0.0, 0.0), Sample::new(10.0, 10.0)];
        let trig = EventSpec::FracMax {
            var: &s,
            frac: 0.1,
            window: Window::FULL,
            edge: Edge::Rising,
            occurrence: Occurrence::Nth(1),
        };
        let targ = EventSpec::FracMax {
            var: &s,
            frac: 0.9,
            window: Window::FULL,
            edge: Edge::Rising,
            occurrence: Occurrence::Nth(1),
        };
        let dt = trig_targ(&trig, &targ).unwrap();
        assert!((dt - 8.0).abs() < 1e-9, "got {dt}");
    }

    /// A crossing between two *different* variables, sampled on different time grids: var1 is a
    /// ramp 0->10 over [0,10], var2 is constant 5. Their difference crosses zero at t=5 (where
    /// var1=5), independent of var2's own sample grid.
    #[test]
    fn crossing_between_two_variables_on_different_grids() {
        let var1 = [Sample::new(0.0, 0.0), Sample::new(10.0, 10.0)];
        let var2 = [
            Sample::new(0.0, 5.0),
            Sample::new(3.0, 5.0),
            Sample::new(10.0, 5.0),
        ];
        let t = when(
            &var1,
            Threshold::Variable(&var2),
            &Window::FULL,
            Edge::Rising,
            Occurrence::Nth(1),
        )
        .unwrap();
        assert!((t - 5.0).abs() < 1e-9, "got {t}");
    }
}
