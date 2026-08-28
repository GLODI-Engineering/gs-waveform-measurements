//! Cycle-counting measurements: `FREQ`, `ON_TIME`, `OFF_TIME`.
//!
//! These follow a two-level (Schmitt-trigger) detector, matching Xyce's description: a signal
//! is considered to enter the "on" state when it rises through `on`, and the "off" state when
//! it falls through `off` (`on` is normally >= `off`; the gap between them is the hysteresis
//! band, which for a clean square/pulse wave can be zero). A "cycle" runs from one on-crossing
//! to the next. `FREQ` is `(number of complete cycles) / (time spanned by those cycles)`.
//! `ON_TIME`/`OFF_TIME` are the total time spent above `on` / below `off` within the window,
//! each normalized (divided) by the number of complete cycles, per Xyce's own definition ("...
//! normalized by the number of cycles of the waveform").
//!
//! Every crossing time used here is found the same way as [`crate::crossing`]: linear
//! interpolation between the two adjacent real samples that straddle the threshold.

use crate::sample::{clip, MeasureError, Result, Sample, Window};

#[derive(Clone, Copy, PartialEq)]
enum State {
    On,
    Off,
    Unknown,
}

struct CycleWalk {
    on_crossings: Vec<f64>,
    /// (on_time_start, on_time_end) intervals where the signal was latched "on".
    on_intervals: Vec<(f64, f64)>,
    /// (off_time_start, off_time_end) intervals where the signal was latched "off".
    off_intervals: Vec<(f64, f64)>,
}

fn walk(samples: &[Sample], on: f64, off: f64) -> CycleWalk {
    let mut on_crossings = Vec::new();
    let mut on_intervals = Vec::new();
    let mut off_intervals = Vec::new();

    let mut state = if samples[0].v >= on {
        State::On
    } else if samples[0].v <= off {
        State::Off
    } else {
        State::Unknown
    };
    let mut state_start = samples[0].t;
    if state == State::On {
        on_crossings.push(samples[0].t);
    }

    for w in samples.windows(2) {
        let (a, b) = (w[0], w[1]);
        if state != State::On && b.v >= on && a.v < on {
            // Rising crossing of `on`.
            let t = a.t + (on - a.v) * (b.t - a.t) / (b.v - a.v);
            if state == State::Off {
                off_intervals.push((state_start, t));
            }
            on_crossings.push(t);
            state = State::On;
            state_start = t;
        } else if state != State::Off && b.v <= off && a.v > off {
            // Falling crossing of `off`.
            let t = a.t + (off - a.v) * (b.t - a.t) / (b.v - a.v);
            if state == State::On {
                on_intervals.push((state_start, t));
            }
            state = State::Off;
            state_start = t;
        }
    }
    let last_t = samples[samples.len() - 1].t;
    match state {
        State::On => on_intervals.push((state_start, last_t)),
        State::Off => off_intervals.push((state_start, last_t)),
        State::Unknown => {}
    }

    CycleWalk {
        on_crossings,
        on_intervals,
        off_intervals,
    }
}

fn cycle_count_and_span(walk: &CycleWalk) -> Result<(usize, f64)> {
    if walk.on_crossings.len() < 2 {
        return Err(MeasureError::NotFound);
    }
    let n = walk.on_crossings.len() - 1;
    let span = walk.on_crossings[walk.on_crossings.len() - 1] - walk.on_crossings[0];
    Ok((n, span))
}

/// `FREQ`: an estimate of the signal's frequency by counting complete on/off cycles between the
/// `on` and `off` thresholds within `window`.
///
/// $$\mathrm{FREQ} = \frac{\text{number of complete cycles}}{t_{\text{last on-crossing}} -
/// t_{\text{first on-crossing}}}$$
pub fn freq(samples: &[Sample], window: &Window, on: f64, off: f64) -> Result<f64> {
    let clipped = clip(samples, window)?;
    let w = walk(&clipped, on, off);
    let (n, span) = cycle_count_and_span(&w)?;
    if span <= 0.0 {
        return Err(MeasureError::NotFound);
    }
    Ok(n as f64 / span)
}

/// Sums each interval's overlap with `[lo, hi]`, discarding any portion outside it. `ON_TIME`/
/// `OFF_TIME` clip to the span between the first and last `on`-crossing (rather than counting
/// the whole window) so that a partial pulse straddling the very start or end of the data —
/// which was never bounded by a full cycle on both sides — doesn't distort the per-cycle
/// average.
fn sum_intervals_clipped(intervals: &[(f64, f64)], lo: f64, hi: f64) -> f64 {
    intervals
        .iter()
        .map(|&(a, b)| (b.min(hi) - a.max(lo)).max(0.0))
        .sum()
}

/// `ON_TIME`: total time `samples` spends at or above `on`, restricted to the span covered by
/// complete cycles, normalized by the number of complete cycles (Xyce's own definition: "the
/// time that `<variable>` is above ON ... normalized by the number of cycles of the waveform").
pub fn on_time(samples: &[Sample], window: &Window, on: f64, off: f64) -> Result<f64> {
    let clipped = clip(samples, window)?;
    let w = walk(&clipped, on, off);
    let (n, span_lo, span_hi) = cycle_count_and_bounds(&w)?;
    Ok(sum_intervals_clipped(&w.on_intervals, span_lo, span_hi) / n as f64)
}

/// `OFF_TIME`: total time `samples` spends at or below `off`, restricted to the span covered by
/// complete cycles, normalized by the number of complete cycles (Xyce's own definition).
pub fn off_time(samples: &[Sample], window: &Window, on: f64, off: f64) -> Result<f64> {
    let clipped = clip(samples, window)?;
    let w = walk(&clipped, on, off);
    let (n, span_lo, span_hi) = cycle_count_and_bounds(&w)?;
    Ok(sum_intervals_clipped(&w.off_intervals, span_lo, span_hi) / n as f64)
}

fn cycle_count_and_bounds(walk: &CycleWalk) -> Result<(usize, f64, f64)> {
    let (n, _) = cycle_count_and_span(walk)?;
    Ok((
        n,
        walk.on_crossings[0],
        walk.on_crossings[walk.on_crossings.len() - 1],
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A clean square wave, period 1: high (=1) on [0,0.5), low (=0) on [0.5,1), repeated 4
    /// times, built from exact PWL segments with near-vertical transitions at dt=1e-9 so the
    /// on/off crossing times are (by hand) essentially exactly at the half-period boundaries.
    /// on=0.5, off=0.5 (no hysteresis, sharp transition threshold at the midpoint).
    fn square_wave(periods: usize) -> Vec<Sample> {
        let mut s = Vec::new();
        let eps = 1e-9;
        let mut t = 0.0;
        s.push(Sample::new(t, 0.0));
        for _ in 0..periods {
            s.push(Sample::new(t + eps, 1.0));
            t += 0.5;
            s.push(Sample::new(t - eps, 1.0));
            s.push(Sample::new(t + eps, 0.0));
            t += 0.5;
            s.push(Sample::new(t - eps, 0.0));
        }
        s.push(Sample::new(t, 0.0));
        s
    }

    /// 4 periods of a 1 Hz square wave: FREQ must be ~1.0 Hz by hand (period=1s).
    #[test]
    fn freq_of_square_wave() {
        let s = square_wave(4);
        let f = freq(&s, &Window::FULL, 0.5, 0.5).unwrap();
        assert!((f - 1.0).abs() < 1e-4, "got {f}");
    }

    /// Same square wave: ON_TIME and OFF_TIME must each be ~0.5s per cycle, by hand (50% duty).
    #[test]
    fn on_off_time_of_square_wave() {
        let s = square_wave(4);
        let on_t = on_time(&s, &Window::FULL, 0.5, 0.5).unwrap();
        let off_t = off_time(&s, &Window::FULL, 0.5, 0.5).unwrap();
        assert!((on_t - 0.5).abs() < 1e-4, "on_time got {on_t}");
        assert!((off_t - 0.5).abs() < 1e-4, "off_time got {off_t}");
    }

    /// A 25%-duty pulse train (high for 0.25s, low for 0.75s, period 1s), 3 periods: ON_TIME
    /// should be ~0.25 and OFF_TIME ~0.75 per cycle, by hand.
    #[test]
    fn asymmetric_duty_cycle() {
        let mut s = Vec::new();
        let eps = 1e-9;
        let mut t = 0.0;
        s.push(Sample::new(t, 0.0));
        for _ in 0..3 {
            s.push(Sample::new(t + eps, 1.0));
            t += 0.25;
            s.push(Sample::new(t - eps, 1.0));
            s.push(Sample::new(t + eps, 0.0));
            t += 0.75;
            s.push(Sample::new(t - eps, 0.0));
        }
        s.push(Sample::new(t, 0.0));
        let on_t = on_time(&s, &Window::FULL, 0.5, 0.5).unwrap();
        let off_t = off_time(&s, &Window::FULL, 0.5, 0.5).unwrap();
        assert!((on_t - 0.25).abs() < 1e-4, "on_time got {on_t}");
        assert!((off_t - 0.75).abs() < 1e-4, "off_time got {off_t}");
    }
}
