# gs-waveform-measurements

`gs-waveform-measurements` is a post-simulation waveform measurement library — max/min/RMS/
integral, threshold crossings, TRIG-TARG, bounded-harmonic Fourier analysis, and error norms —
over a generic `(time, value)` time series, in the style of ngspice's `.measure` and Xyce's
`.MEASURE`, computed exactly on the real samples rather than by resampling to a uniform grid.

[![crates.io](https://img.shields.io/crates/v/gs-waveform-measurements.svg)](https://crates.io/crates/gs-waveform-measurements)
[![docs.rs](https://img.shields.io/docsrs/gs-waveform-measurements)](https://docs.rs/gs-waveform-measurements)

```toml
[dependencies]
gs-waveform-measurements = "0.1"
```

The API reference is at
[docs.rs/gs-waveform-measurements](https://docs.rs/gs-waveform-measurements).

## Why this exists

A variable-timestep transient solver's own output is naturally non-uniformly spaced in time —
the solver only records the solution at the instants it actually stepped to. Measurement code
that is not written with that in mind gets subtly wrong answers: the one wrong thing to
implement is a plain `sum(values) / count` for anything described as an "average" — that
silently misweights variable step sizes. This crate follows the three techniques ngspice's
`.measure`, Xyce's `.MEASURE`, and a block-diagram simulation tool's own documented
Fourier-coefficient method each use to get this right on real (non-uniform) sample data:

1. **Scalar measurements** (`RMS`/`AVG`/`INTEG`/...) are Δt-weighted — trapezoidal integration
   over the real samples, never a naive mean.
2. **Crossing/threshold measurements** (`WHEN`, `FIND ... WHEN`, `TRIG`/`TARG`, ...) are found by
   linear interpolation strictly between the two adjacent real samples that straddle the
   threshold — never by resampling the whole trace.
3. **Bounded-harmonic Fourier analysis** (`FOUR`) is computed by exact analytic integration of
   each linear segment against the harmonic kernel, summed in closed form — not a resampled FFT.

## Design constraint: no circuit knowledge

This crate takes only `Sample` slices as input and returns plain numeric results. It has no
dependency on any circuit or block-diagram solver — it is meant to be just as usable from a
project with no connection to any particular simulator at all.

It deliberately does not read a `.measure`/`.MEASURE` netlist statement, does not know about a
CLI, and does not read a "comparison waveform" out of a file for the `ERROR` measure — wiring a
text-format measurement syntax and file I/O onto these primitives is a separate, later
integration layer left to a consumer. Two measure types from the ngspice/Xyce catalogue are
intentionally not implemented at all here: `EQN` (evaluating an expression over *other
measurements'* results is an orchestration-layer concern, not a per-signal waveform primitive)
and the `FILE=` half of `ERROR` (parsing a reference-waveform file format). The norm computation
itself for `ERROR` — comparing two given time series — *is* implemented (`error_norm`).

## Basic use

```rust
use gs_waveform_measurements::{integ, rms, Sample, Window};

let samples = [
    Sample::new(0.0, 0.0),
    Sample::new(0.5, 10.0),
    Sample::new(1.5, -4.0),
    Sample::new(3.0, 2.0),
];

// INTEG: trapezoidal integral over the real (non-uniformly spaced) samples.
let area = integ(&samples, &Window::FULL).expect("series has at least two samples");

// RMS: exact for piecewise-linear data (each segment's v^2 integrated in closed form,
// not approximated by trapezoidally integrating the squared sample values).
let effective = rms(&samples, &Window::FULL).expect("series has at least two samples");
```

`Window` mirrors ngspice/Xyce's `FROM=`/`TO=`/`TD=` qualifiers — pass `Window::new(from, to)` to
restrict a measurement to part of the series; `Window::FULL` covers the whole thing.

## The six modules

| module | what it measures | ngspice/Xyce equivalent |
| --- | --- | --- |
| `scalar` | `max`, `min`, `peak_to_peak`, `integ`, `avg`, `rms` over a window | `MAX`/`MIN`/`PP`/`INTEG`/`AVG`/`RMS` |
| `crossing` | `find_at`, `when`, `find_when`, `deriv_when`, `trig_targ` | `FIND ... AT=`, `WHEN`, `TRIG ... TARG` |
| `cycle` | `freq`, `on_time`, `off_time` from a two-level signal | `FREQ`, `ON_TIME`, `OFF_TIME` |
| `fourier` | `four` — bounded-harmonic magnitudes, phases and THD | `FOUR` with `NUMFREQ` |
| `error_metrics` | `err1`, `err2`, `error_norm` between two series | `ERROR` with `COMP_FUNCTION=` |
| `sample` | `Sample`, `Window`, `interpolate`, `clip`, `check_monotonic` | `FROM=`/`TO=`/`TD=` |

The three that are easiest to get wrong are worth a worked example each. All three are compiled
and checked against hand-derived values in
[`tests/readme_examples.rs`](tests/readme_examples.rs), so they cannot drift from the API.

### `trig_targ`: a rise time between two fractions of the signal's own peak

```rust
use gs_waveform_measurements::{trig_targ, Edge, EventSpec, Occurrence, Sample, Window};

// A step, sampled at the instants a variable-step solver happened to stop at.
let step = [
    Sample::new(0.0, 0.0),
    Sample::new(1.0e-6, 0.0),
    Sample::new(1.4e-6, 0.8),
    Sample::new(2.0e-6, 1.0),
];
let at_fraction = |frac| EventSpec::FracMax {
    var: &step,
    frac,
    window: Window::FULL,
    edge: Edge::Rising,
    occurrence: Occurrence::Nth(1),
};

// 650 ns. Neither the 10 % nor the 90 % level falls on a sample: each is interpolated
// strictly between the two real samples that straddle it.
let t_rise = trig_targ(&at_fraction(0.1), &at_fraction(0.9))?;
# Ok::<(), gs_waveform_measurements::MeasureError>(())
```

`FracMax` takes the level from `max(var)` over the window, which is what you want when the
signal's final value is not known in advance. Use `EventSpec::Crossing` for an absolute level and
`EventSpec::At` for a time you already know.

### `four`: harmonics and THD, exactly, on piecewise-linear data

```rust
use gs_waveform_measurements::{four, Sample, Window};

let period = 0.02; // 50 Hz
let triangle = [
    Sample::new(0.0, 0.0),
    Sample::new(period / 4.0, 1.0),
    Sample::new(3.0 * period / 4.0, -1.0),
    Sample::new(period, 0.0),
];

// DC plus orders 1..4 — ngspice/Xyce's `NUMFREQ` convention, DC included.
let spectrum = four(&triangle, &Window::new(0.0, period), 50.0, 5)?;

// A unit triangle's harmonics are 8/(pi^2 n^2) for odd n and zero for even n, and a triangle
// is exactly piecewise-linear, so these come out exact rather than approximate:
// harmonics[1].magnitude == 8/pi^2, harmonics[3].magnitude == 8/(9 pi^2), harmonics[2] == 0.
// THD counts only the harmonics inside the bound, so here it is exactly 100/9 %.
assert!((spectrum.thd_percent - 100.0 / 9.0).abs() < 1e-9);
# Ok::<(), gs_waveform_measurements::MeasureError>(())
```

`window`'s duration is the integration period: choosing a window that spans a whole number of
fundamental periods is the caller's job, exactly as in ngspice and Xyce.

### `error_norm`: how far a run is from a reference waveform

```rust
use gs_waveform_measurements::{error_norm, Norm, Sample};

let measured = [Sample::new(0.0, 0.0), Sample::new(1.0, 1.0), Sample::new(2.0, 2.0)];
let reference = [Sample::new(0.0, 0.0), Sample::new(1.0, 1.5), Sample::new(2.0, 2.0)];

// Evaluated at `measured`'s own sample times, over the overlap of both domains, with
// `reference` interpolated onto that grid: sqrt(0.25/3) = 0.2887.
let l2 = error_norm(&measured, &reference, Norm::L2)?;
# Ok::<(), gs_waveform_measurements::MeasureError>(())
```

`Norm::L1`, `Norm::L2` (Xyce's default) and `Norm::InfNorm` are the three Xyce `ERROR` norms.
Reading the reference out of a file is deliberately not this crate's job.

## Running the tests

```bash
cargo test --all-targets
```

Every measurement's test suite is verified against hand-derived expected values (closed-form
triangle areas, analytic sine RMS, textbook Fourier coefficients, ...), not merely internal
self-consistency — see `AGENTS.md` for this project's verification discipline.

## Development

Quality gates:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --all-targets
cargo doc --no-deps --lib
```

`cargo doc` should stay clean — the crate is built with `#![warn(missing_docs)]`.

### Who depends on this

Two consumers, and both are why the crate is published rather than vendored:

- **[general-simulator](https://github.com/GLODI-Engineering/general-simulator)** — the engine's
  `kind=measure` netlist cards are evaluated with these primitives, and its windowed rules are
  being built on them. Its CI, and `general-mna`'s, check this repository out as a sibling for
  exactly that reason.
- **[gsim-tools](https://github.com/GLODI-Engineering/gsim-tools)** — `gsim_tools.measure` is a
  Python port of this crate, cross-validated against it measurement by measurement in that
  repository's `tests/test_measure_parity.py`. A measurement that disagrees between the two is a
  failing test there.

During development on a checkout that has the siblings beside it, a path dependency is the
convention (the same one `general-mna`'s README documents); a build from crates.io resolves the
published version:

```toml
gs-waveform-measurements = { path = "../gs-waveform-measurements", version = "0.1.0" }
```

### History

This crate was originally developed inside the `general-simulator` workspace and extracted into
its own repository so it can carry an independent release cadence and be reused outside that one
project. See `general-simulator`'s own `docs/journal/` for the crate's original design rationale
and non-uniform-timestep verification strategy, and this repository's own `docs/journal/` for the
extraction itself.

Recurring workflows live in [`.claude/skills`](.claude/skills), and the dated continuity log
lives in [`docs/journal`](docs/journal).
