# Agent guidelines

`gs-waveform-measurements` is a small Rust library of post-simulation waveform measurements
(max/min/RMS/integral, threshold crossings, TRIG-TARG, bounded-harmonic Fourier analysis, error
norms) over a generic, non-uniformly-spaced `(time, value)` series. It has no dependency on any
circuit or block-diagram solver.

## Authority order

1. Mechanical gates: formatters, Clippy, tests, `cargo doc`.
2. Skills in `.claude/skills/<name>/SKILL.md` for recurring workflows.
3. This file.
4. `docs/` (journal, gotchas) for durable rationale.

Never bypass a failing gate. Fix the cause.

## Project boundaries

- No dependency on any circuit, block-diagram, or netlist type — this crate takes only `Sample`
  slices and returns plain numeric results, so it stays usable from a project with no connection
  to any particular simulator.
- Every measurement operates on the real, possibly non-uniformly-spaced samples it is given.
  Never resample to a uniform grid; never implement an "average" as a plain
  `sum(values) / count` — that silently misweights variable step sizes. See `src/lib.rs`'s
  crate-level doc comment for the three techniques (Δt-weighted trapezoidal/closed-form
  integration, straddling-sample linear interpolation, exact per-segment analytic Fourier
  integration) every measurement in this crate follows.
- Describe a technique adapted from a named external tool generically (e.g. "a block-diagram
  simulation tool") rather than naming it — see `scripts/style-guard.sh` if present.

## Layout

```text
src/          measurement modules: sample, scalar, crossing, cycle, fourier, error_metrics
tests/        Rust integration tests (unit tests also live inline per module)
docs/         journal and gotchas
```

## Development cycle

The gates below say how code is written. This says how the work is recorded, which matters for
the same reason: someone reads it later, and what they can reconstruct is what was written down.

1. **A large instruction becomes a milestone**, named for the outcome rather than the task
   ("Rules and diagnostics", not "implement the rules spec").
2. **One issue per deliverable phase or decision**, each stating the problem rather than the
   patch. The instruction as a whole — its reasoning, what it leaves open — is an issue in that
   milestone too, and the rest reference it with `Part of #N`.
3. **A branch and a pull request per issue**, the description ending `Closes #N`, assigned to the
   milestone. Nothing lands on the default branch directly.
4. **The merge closes the issue.** Closing it by hand loses the link between the issue and the
   commit that answered it, which is the whole point of the exercise.
5. **Anything found but not done becomes its own issue**, not a sentence in a pull request body.
   A merged pull request's prose is invisible; an open issue is not. A deliberate omission with
   a reason is worth more than a half-finished corner of a change.
6. **Gates run before pushing**, never after. A red branch spends someone's CI and tells the
   reader the author did not look.

Finish a repository's open pull request before starting new work in it, so review has one thing
to hold at a time.

## Required gates

Before committing, use the `commit` skill (if present) and run:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --all-targets
cargo doc --no-deps --lib
```

`cargo doc` must stay clean — the crate is built with `#![warn(missing_docs)]`.

## Verification discipline

Numerical agreement with the code's own internal self-consistency is not proof. Every
measurement's tests are checked against a hand-derived expected value — a closed-form triangle
area for `INTEG`, the analytic `amplitude / sqrt(2)` for a sine's RMS, textbook Fourier
coefficients, a known crossing time from linear interpolation by hand — not just "the function
returns a plausible-looking number" or agreement between two code paths in this crate. Any new
measurement should follow the same discipline: derive the expected value independently before
writing the assertion.

## Skills

| Task | Skill |
|---|---|
| Gate, review, and commit changes | `.claude/skills/commit/SKILL.md` |

## Journal and gotchas

Read the newest top entry in `docs/journal/` when continuing prior work. Include one dated
journal entry with every commit, using a time verified with `date '+%Y-%m-%d %H:%M'`. Keep
entries newest-first and append history rather than rewriting it.

Use `docs/gotchas/` (see `_TEMPLATE.md`) for a reproducible trap another contributor would
otherwise have to rediscover.

## Commit rules

- Conventional Commits: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `build`,
  `ci`, `chore`, or `revert`.
- Review the complete diff and stage deliberately.
- Never use `--no-verify`, `SKIP=`, or similar bypasses.
