//! Scalar reductions over a window: `MAX`/`MIN`/`MAX_AT`/`MIN_AT`, `PP`, `AVG`, `RMS`, `INTEG`.
//!
//! None of these resample the input. `INTEG` and `AVG` use the trapezoidal rule, which is
//! *exact* (not an approximation) for a signal that really is piecewise-linear between the
//! given samples — this is the same method Xyce's own reference guide describes for `INTEG`
//! ("second order numerical integration"). `RMS` goes one step further: rather than
//! trapezoidally integrating the *squared* sample values (which would only approximate
//! $\int v^2\,dt$, since $v^2$ is quadratic, not linear, within a segment), it integrates each
//! linear segment's square in closed form, so it is exact for PWL data too. Every one of these
//! is $\Delta t$-weighted — none of them is a plain `sum(values) / count`, which silently
//! misweights variable step sizes.

use crate::sample::{clip, Result, Sample, Window};

/// The value and time of an extremum (`MAX`/`MIN`), covering both ngspice's `MAX_AT`/`MIN_AT`
/// and Xyce's `MAX ... OUTPUT=TIME|VALUE` in a single result: read `.value` for the ngspice-style
/// `MAX`/`MIN` result, `.time` for `MAX_AT`/`MIN_AT` or Xyce's `OUTPUT=TIME`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extremum {
    /// The extreme value itself.
    pub value: f64,
    /// The time at which that value occurs. If the extremum is attained on a flat plateau, this
    /// is the time of the *first* sample (scanning forward) at which it occurs.
    pub time: f64,
}

/// `MAX` / `MAX_AT`: the maximum value of `samples` within `window`, and the time it occurs at.
///
/// $$\mathrm{MAX} = \max_i v_i \quad \text{over } t_i \in \text{window (with interpolated
/// boundary samples)}$$
pub fn max(samples: &[Sample], window: &Window) -> Result<Extremum> {
    let clipped = clip(samples, window)?;
    let best = clipped
        .into_iter()
        .fold(None::<Extremum>, |acc, s| match acc {
            Some(e) if e.value >= s.v => Some(e),
            _ => Some(Extremum {
                value: s.v,
                time: s.t,
            }),
        });
    Ok(best.expect("clip() guarantees at least two samples"))
}

/// `MIN` / `MIN_AT`: the minimum value of `samples` within `window`, and the time it occurs at.
pub fn min(samples: &[Sample], window: &Window) -> Result<Extremum> {
    let clipped = clip(samples, window)?;
    let best = clipped
        .into_iter()
        .fold(None::<Extremum>, |acc, s| match acc {
            Some(e) if e.value <= s.v => Some(e),
            _ => Some(Extremum {
                value: s.v,
                time: s.t,
            }),
        });
    Ok(best.expect("clip() guarantees at least two samples"))
}

/// `PP`: peak-to-peak amplitude, `max(v) - min(v)`, within `window`.
pub fn peak_to_peak(samples: &[Sample], window: &Window) -> Result<f64> {
    let hi = max(samples, window)?;
    let lo = min(samples, window)?;
    Ok(hi.value - lo.value)
}

/// `INTEG`: the definite integral $\int_{t_0}^{t_1} v\,dt$ over `window`, via the trapezoidal
/// rule applied directly to the (non-uniformly spaced) real samples — exact for PWL data,
/// correct under variable time steps because each trapezoid uses that segment's own $\Delta t$.
///
/// $$\mathrm{INTEG} = \sum_i \frac{v_i + v_{i+1}}{2}\,(t_{i+1} - t_i)$$
pub fn integ(samples: &[Sample], window: &Window) -> Result<f64> {
    let clipped = clip(samples, window)?;
    Ok(trapezoid(&clipped))
}

fn trapezoid(samples: &[Sample]) -> f64 {
    samples
        .windows(2)
        .map(|w| 0.5 * (w[0].v + w[1].v) * (w[1].t - w[0].t))
        .sum()
}

/// `AVG`: the time ($\Delta t$-)weighted mean of `samples` over `window`, i.e. `INTEG / duration`
/// — never a plain arithmetic mean of sample values, which would misweight irregular spacing.
///
/// $$\mathrm{AVG} = \frac{1}{t_1 - t_0}\int_{t_0}^{t_1} v\,dt$$
pub fn avg(samples: &[Sample], window: &Window) -> Result<f64> {
    let clipped = clip(samples, window)?;
    let duration = clipped[clipped.len() - 1].t - clipped[0].t;
    Ok(trapezoid(&clipped) / duration)
}

/// `RMS`: the $\Delta t$-weighted root-mean-square of `samples` over `window`.
///
/// Xyce's reference guide defines RMS as "the square root of the area under the `<variable>`
/// curve, divided by the period of interest" (understood as the area under $v^2$). Since $v$
/// is linear within each segment
/// `[t_i, t_{i+1}]` with $v(u) = v_i + m\,u$ for $u \in [0, \Delta t]$, $v^2$ is a quadratic in
/// $u$ and its exact integral over the segment is:
///
/// $$\int_0^{\Delta t} v(u)^2\,du = \frac{\Delta t}{3}\left(v_i^2 + v_i v_{i+1} +
/// v_{i+1}^2\right)$$
///
/// (expand $v(u)^2 = v_i^2 + 2 v_i m u + m^2 u^2$ and integrate term-by-term; the constant
/// $2v_i m \Delta t^2/2 + m^2\Delta t^3/3$ terms simplify to this symmetric form once
/// $m\,\Delta t = v_{i+1}-v_i$ is substituted). Summing this over every segment and dividing by
/// the total duration gives an RMS that is *exact* for PWL data — not merely an approximation
/// from trapezoidally integrating the squared sample values, which would use only a linear
/// (not quadratic) fit to $v^2$ within each segment.
pub fn rms(samples: &[Sample], window: &Window) -> Result<f64> {
    let clipped = clip(samples, window)?;
    let duration = clipped[clipped.len() - 1].t - clipped[0].t;
    let integral_v2: f64 = clipped
        .windows(2)
        .map(|w| {
            let dt = w[1].t - w[0].t;
            let (v0, v1) = (w[0].v, w[1].v);
            dt / 3.0 * (v0 * v0 + v0 * v1 + v1 * v1)
        })
        .sum();
    Ok((integral_v2 / duration).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A single ramp from (0,0) to (1,1): INTEG = area of a right triangle = 1/2, by hand.
    #[test]
    fn integ_ramp_matches_triangle_area() {
        let s = [Sample::new(0.0, 0.0), Sample::new(1.0, 1.0)];
        assert!((integ(&s, &Window::FULL).unwrap() - 0.5).abs() < 1e-12);
    }

    /// Non-uniform Δt: three points 0,0 -> 1,10 -> 4,10 (dt=1 then dt=3). Trapezoid areas by
    /// hand: seg1 = (0+10)/2*1 = 5; seg2 = (10+10)/2*3 = 30. Total = 35. AVG = 35/4 = 8.75.
    /// A naive sum(values)/count = (0+10+10)/3 = 6.667 would be wrong — this proves the
    /// implementation is Δt-weighted, not a plain mean.
    #[test]
    fn avg_is_dt_weighted_not_plain_mean() {
        let s = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 10.0),
            Sample::new(4.0, 10.0),
        ];
        let a = avg(&s, &Window::FULL).unwrap();
        assert!((a - 8.75).abs() < 1e-12, "got {a}");
        assert!(
            (a - 6.666_666_666_666_667).abs() > 0.1,
            "must not equal naive mean"
        );
    }

    /// A pure sine wave over exactly one period: RMS = amplitude / sqrt(2), analytically. Use a
    /// deliberately non-uniform sampling (irregular Δt) to prove the Δt-weighting, not just a
    /// convenient uniform grid.
    #[test]
    fn rms_of_sine_matches_amplitude_over_sqrt2_under_nonuniform_sampling() {
        let amplitude = 3.0;
        let freq = 1.0; // Hz
        let period = 1.0 / freq;
        // Irregular time steps: not evenly spaced, but dense enough that the PWL approximation
        // of a smooth sine is accurate to the asserted tolerance.
        let mut ts = vec![0.0];
        let mut t = 0.0;
        let mut i = 0usize;
        while t < period {
            let dt = 0.0005 + 0.0004 * ((i as f64) * 0.7).sin().abs();
            t += dt;
            if t > period {
                t = period;
            }
            ts.push(t);
            i += 1;
        }
        let samples: Vec<Sample> = ts
            .iter()
            .map(|&t| Sample::new(t, amplitude * (2.0 * std::f64::consts::PI * freq * t).sin()))
            .collect();
        let r = rms(&samples, &Window::FULL).unwrap();
        let expected = amplitude / std::f64::consts::SQRT_2;
        assert!((r - expected).abs() < 1e-3, "got {r}, expected {expected}");
    }

    /// A signal constant at value `c` over duration `d`: RMS = |c| exactly, regardless of
    /// (non-uniform) sample placement, since v^2 is exactly c^2 pointwise. This exercises the
    /// exact-quadrature RMS formula on a case with a trivially known closed-form answer.
    #[test]
    fn rms_of_constant_is_exact_regardless_of_spacing() {
        let s = [
            Sample::new(0.0, 5.0),
            Sample::new(0.1, 5.0),
            Sample::new(0.15, 5.0),
            Sample::new(1.0, 5.0),
        ];
        let r = rms(&s, &Window::FULL).unwrap();
        assert!((r - 5.0).abs() < 1e-12, "got {r}");
    }

    #[test]
    fn max_min_pp_and_at_times() {
        let s = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 10.0),
            Sample::new(2.0, -4.0),
            Sample::new(3.0, 2.0),
        ];
        let mx = max(&s, &Window::FULL).unwrap();
        assert_eq!(
            mx,
            Extremum {
                value: 10.0,
                time: 1.0
            }
        );
        let mn = min(&s, &Window::FULL).unwrap();
        assert_eq!(
            mn,
            Extremum {
                value: -4.0,
                time: 2.0
            }
        );
        assert!((peak_to_peak(&s, &Window::FULL).unwrap() - 14.0).abs() < 1e-12);
    }

    /// Window [1.5, 3.0] clips the series, inserting an interpolated boundary sample at t=1.5:
    /// linear interpolation between (1,10) and (2,-4) gives v = 10 + 0.5*(-4-10) = 3.0, by hand
    /// — which is itself the new maximum over the clipped window (larger than the interior
    /// sample -4.0 and the endpoint 2.0), proving max() respects the interpolated boundary.
    #[test]
    fn window_restricts_extrema() {
        let s = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 10.0),
            Sample::new(2.0, -4.0),
            Sample::new(3.0, 2.0),
        ];
        let w = Window::new(1.5, 3.0);
        let mx = max(&s, &w).unwrap();
        assert!((mx.value - 3.0).abs() < 1e-12, "got {}", mx.value);
        assert!((mx.time - 1.5).abs() < 1e-12, "got {}", mx.time);
    }

    #[test]
    fn empty_series_errors() {
        let s: [Sample; 0] = [];
        assert!(max(&s, &Window::FULL).is_err());
    }
}
