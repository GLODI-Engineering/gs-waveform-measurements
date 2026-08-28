# gs-waveform-measurements

`gs-waveform-measurements` is a post-simulation waveform measurement library — max/min/RMS/
integral, threshold crossings, TRIG-TARG, bounded-harmonic Fourier analysis, and error norms —
over a generic `(time, value)` time series, in the style of ngspice's `.measure` and Xyce's
`.MEASURE`, computed exactly on the real samples rather than by resampling to a uniform grid.

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

See the crate's own `cargo doc` output for the full API: scalar reductions (`scalar` module),
crossing/threshold/cycle measurements (`crossing`, `cycle`), bounded-harmonic Fourier analysis
(`fourier`), and error norms between two series (`error_metrics`).

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

### History

This crate was originally developed inside the `general-simulator` workspace and extracted into
its own repository so it can carry an independent release cadence and be reused outside that one
project. See `general-simulator`'s own `docs/journal/` for the crate's original design rationale
and non-uniform-timestep verification strategy, and this repository's own `docs/journal/` for the
extraction itself. Once a consumer in a sibling repository needs it again, it is expected to be
added back as a sibling-path dependency during development, e.g.:

```toml
gs-waveform-measurements = { path = "../gs-waveform-measurements" }
```

(the same sibling-path-during-development convention `general-mna`'s own README documents),
switching to a published `version` once this crate has one.

Recurring workflows live in [`.claude/skills`](.claude/skills), and the dated continuity log
lives in [`docs/journal`](docs/journal).
