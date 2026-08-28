#![warn(missing_docs)]
//! Post-simulation waveform measurements — an ngspice `.measure` / Xyce `.MEASURE`-style
//! measurement library over a generic `(time, value)` time series.
//!
//! # Why this exists
//!
//! A variable-timestep transient solver's own output is naturally non-uniformly spaced in time
//! — the solver only records the solution at the instants it actually stepped to. Measurement
//! code that is not written with that in mind gets subtly wrong answers: the one wrong thing to
//! implement is a plain `sum(values) / count` for anything described as an "average" — that
//! silently misweights variable step sizes. This crate follows the three techniques ngspice's
//! `.measure`, Xyce's `.MEASURE`, and a block-diagram simulation tool's own documented
//! Fourier-coefficient method each use to get this right on real (non-uniform) sample data:
//!
//! 1. **Scalar measurements** (`RMS`/`AVG`/`INTEG`/...) are $\Delta t$-weighted — trapezoidal
//!    integration over the real samples, never a naive mean. See the `scalar` module.
//! 2. **Crossing/threshold measurements** (`WHEN`, `FIND ... WHEN`, `TRIG`/`TARG`, ...) are
//!    found by linear interpolation strictly between the two adjacent real samples that
//!    straddle the threshold — never by resampling the whole trace. See the `crossing` and
//!    `cycle` modules.
//! 3. **Bounded-harmonic Fourier analysis** (`FOUR`) is computed by exact analytic integration
//!    of each linear segment against the harmonic kernel, summed in closed form — not a
//!    resampled FFT. See the `fourier` module.
//!
//! # Design constraint: no circuit knowledge
//!
//! This crate takes only [`Sample`] slices as input and returns plain numeric results. It has
//! **no dependency on this workspace's own circuit or block-diagram types** (`general-mna`-style
//! solvers, `dae-runtime`, `continuous-blocks`, or anything else in this repository) — the same
//! isolation this repository's own `lcp-solver` crate is built and trusted under before any
//! circuit-facing code depends on it. It is meant to be just as usable from a project with no
//! connection to this simulator at all.
//!
//! **What this crate deliberately does *not* do** (see this repository's journal for the full
//! rationale): it does not read a `.measure`/`.MEASURE` netlist statement, does not know about a
//! CLI, and does not read a "comparison waveform" out of a file for the `ERROR` measure —
//! wiring a text-format measurement syntax and file I/O onto these primitives is a separate,
//! later integration layer. Two measure types from the ngspice/Xyce catalogue are intentionally
//! *not implemented at all* here (not merely deferred): `EQN` (evaluating an expression over
//! *other measurements'* results is an orchestration-layer concern, not a per-signal waveform
//! primitive) and the `FILE=` half of `ERROR` (parsing a reference-waveform file format). The
//! norm computation itself for `ERROR` — comparing two given time series — *is* implemented,
//! see [`error_metrics::error_norm`].

mod crossing;
mod cycle;
mod error_metrics;
mod fourier;
mod sample;
mod scalar;

pub use crossing::{
    deriv_at, deriv_when, find_at, find_when, resolve_event, trig_targ, when, Edge, EventSpec,
    Occurrence, Threshold,
};
pub use cycle::{freq, off_time, on_time};
pub use error_metrics::{err1, err2, error_norm, Norm};
pub use fourier::{four, FourierResult, Harmonic};
pub use sample::{check_monotonic, clip, interpolate, MeasureError, Result, Sample, Window};
pub use scalar::{avg, integ, max, min, peak_to_peak, rms, Extremum};
