//! `FOUR`: bounded-harmonic Fourier-coefficient analysis, computed by *exact analytic
//! integration of each piecewise-linear segment* against the harmonic kernel
//! $e^{-j n \omega_0 t}$ — not by resampling the trace onto a uniform grid and running an FFT.
//!
//! This deliberately follows the same technique a block-diagram simulation tool's own
//! documented Fourier-coefficient method uses: because $v(t)$ is linear within each segment
//! `[t_i, t_{i+1}]`, $\int v(t)\,e^{-j n \omega_0 t}\,dt$ has a closed form on that segment, so
//! summing the closed form over every segment gives the exact coefficient for data that really
//! is piecewise-linear between the given (non-uniformly spaced) samples — zero interpolation
//! error, and it costs `O(N * num_harmonics)` rather than the `O(N log N)` of an FFT preceded by
//! a resampling step (and without that resampling step's own error). This is a **bounded-
//! harmonic** measurement, not a broadband FFT: it only ever evaluates a caller-chosen number of
//! harmonics of a caller-chosen fundamental, which is exactly what ngspice's/Xyce's `FOUR`/
//! `.FOUR` compute (a full broadband spectrum is a different, harder feature, out of scope
//! here).

use crate::sample::{clip, MeasureError, Result, Sample, Window};

#[derive(Clone, Copy, Debug)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    const ZERO: Complex = Complex { re: 0.0, im: 0.0 };

    fn add(self, other: Complex) -> Complex {
        Complex {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }

    fn scale(self, k: f64) -> Complex {
        Complex {
            re: self.re * k,
            im: self.im * k,
        }
    }

    fn mul(self, other: Complex) -> Complex {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }

    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
}

/// One harmonic of a [`four`] result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Harmonic {
    /// The harmonic order `n` (0 = DC, 1 = fundamental, 2 = first overtone, ...).
    pub order: usize,
    /// The frequency this harmonic is evaluated at, `n * fundamental_freq`.
    pub frequency: f64,
    /// For `order == 0`: the DC (average) value itself, signed. For `order >= 1`: the peak
    /// amplitude of that harmonic's cosine component, `2 * |C_n|` — always non-negative.
    pub magnitude: f64,
    /// Phase in degrees (matching ngspice/Xyce's convention of reporting `FOUR` phase in
    /// degrees), such that
    /// $v(t) = \mathrm{magnitude}_0 + \sum_{n\ge1} \mathrm{magnitude}_n \cos(n\omega_0 t +
    /// \mathrm{phase}_n)$. Always `0` for the DC term.
    pub phase_degrees: f64,
}

/// The full result of a [`four`] call.
#[derive(Debug, Clone, PartialEq)]
pub struct FourierResult {
    /// One entry per requested harmonic, `orders[0]` is the DC term.
    pub harmonics: Vec<Harmonic>,
    /// Total harmonic distortion, as a percentage:
    /// $$\mathrm{THD} = 100\% \times \frac{\sqrt{\sum_{n \ge 2} A_n^2}}{A_1}$$
    /// (the ratio of the RMS of every harmonic above the fundamental to the fundamental's own
    /// magnitude — the `sqrt(2)` peak-to-RMS factor is common to every term and cancels in the
    /// ratio, so it is equivalent to compute the ratio directly from the peak magnitudes $A_n$).
    /// Requires at least a fundamental and one overtone (`num_harmonics >= 3`, i.e. DC +
    /// fundamental + at least one overtone); otherwise `f64::NAN`.
    pub thd_percent: f64,
}

/// Exact analytic integral of a linear segment against `e^{-j omega t}` for `omega != 0`,
/// derived by direct calculus (see the derivation in this module's own history / journal entry
/// for the closed form) and specialized to `a = j*omega` (purely imaginary) throughout, since
/// that is the only case `FOUR` ever needs — this avoids any general complex-division code path.
fn segment_harmonic_integral(t0: f64, v0: f64, v1: f64, dt: f64, omega: f64) -> Complex {
    let m = (v1 - v0) / dt;
    let (s, c) = (omega * dt).sin_cos();

    // I0 = (1 - e^{-j*omega*dt}) / (j*omega)
    let i0 = Complex {
        re: s / omega,
        im: -(1.0 - c) / omega,
    };
    // I1 = (1 - e^{-j*omega*dt}*(1 + j*omega*dt)) / (j*omega)^2
    let omega2 = omega * omega;
    let i1 = Complex {
        re: (c + omega * dt * s - 1.0) / omega2,
        im: (omega * dt * c - s) / omega2,
    };

    let bracket = i0.scale(v0).add(i1.scale(m));
    // e^{-j*omega*t0}
    let (s0, c0) = (omega * t0).sin_cos();
    let phase = Complex { re: c0, im: -s0 };
    phase.mul(bracket)
}

/// `FOUR`: Fourier-coefficient analysis of `samples` over `window`, at fundamental frequency
/// `fundamental_freq` (Hz), reporting `num_harmonics` terms (DC included, so `num_harmonics = N`
/// yields orders `0..=N-1`, matching ngspice/Xyce's `NUMFREQ` convention of "DC plus
/// `NUMFREQ - 1` harmonics").
///
/// `window`'s duration is used as the integration period `T` (the caller is responsible for
/// choosing a window that corresponds to a meaningful number of fundamental periods, exactly as
/// in ngspice/Xyce — this crate has no notion of "the" period of an arbitrary signal).
pub fn four(
    samples: &[Sample],
    window: &Window,
    fundamental_freq: f64,
    num_harmonics: usize,
) -> Result<FourierResult> {
    if num_harmonics == 0 {
        return Err(MeasureError::InvalidParameter(
            "num_harmonics must be at least 1 (the DC term)",
        ));
    }
    if fundamental_freq <= 0.0 {
        return Err(MeasureError::InvalidParameter(
            "fundamental_freq must be positive",
        ));
    }
    let clipped = clip(samples, window)?;
    let t_duration = clipped[clipped.len() - 1].t - clipped[0].t;

    let mut harmonics = Vec::with_capacity(num_harmonics);
    for n in 0..num_harmonics {
        if n == 0 {
            let dc: f64 = clipped
                .windows(2)
                .map(|w| 0.5 * (w[0].v + w[1].v) * (w[1].t - w[0].t))
                .sum::<f64>()
                / t_duration;
            harmonics.push(Harmonic {
                order: 0,
                frequency: 0.0,
                magnitude: dc,
                phase_degrees: 0.0,
            });
            continue;
        }
        let omega = 2.0 * std::f64::consts::PI * (n as f64) * fundamental_freq;
        let mut total = Complex::ZERO;
        for w in clipped.windows(2) {
            let dt = w[1].t - w[0].t;
            total = total.add(segment_harmonic_integral(w[0].t, w[0].v, w[1].v, dt, omega));
        }
        let c_n = total.scale(1.0 / t_duration);
        harmonics.push(Harmonic {
            order: n,
            frequency: n as f64 * fundamental_freq,
            magnitude: 2.0 * c_n.abs(),
            phase_degrees: c_n.im.atan2(c_n.re).to_degrees(),
        });
    }

    let thd_percent = if num_harmonics >= 3 {
        let fundamental = harmonics[1].magnitude;
        let overtones_sum_sq: f64 = harmonics[2..]
            .iter()
            .map(|h| h.magnitude * h.magnitude)
            .sum();
        100.0 * overtones_sum_sq.sqrt() / fundamental
    } else {
        f64::NAN
    };

    Ok(FourierResult {
        harmonics,
        thd_percent,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// A near-ideal sampled sine wave (2000 points/period, so the PWL discretization error is
    /// negligible) must have its n=1 harmonic magnitude close to the true amplitude and near-
    /// zero THD, by hand (a pure sinusoid has no harmonic content).
    #[test]
    fn pure_sine_has_expected_fundamental_and_low_thd() {
        let amplitude = 2.0;
        let freq = 50.0;
        let period = 1.0 / freq;
        let n_points = 2000;
        let samples: Vec<Sample> = (0..=n_points)
            .map(|i| {
                let t = period * (i as f64) / (n_points as f64);
                Sample::new(t, amplitude * (2.0 * PI * freq * t).sin())
            })
            .collect();
        let result = four(&samples, &Window::FULL, freq, 5).unwrap();
        let fundamental = result.harmonics[1];
        assert!(
            (fundamental.magnitude - amplitude).abs() < 1e-3,
            "got {}",
            fundamental.magnitude
        );
        assert!(result.thd_percent < 0.1, "got {}", result.thd_percent);
        // v(t) = A sin(wt) = A cos(wt - 90deg), so phase should be -90 degrees, by hand.
        assert!(
            (fundamental.phase_degrees - (-90.0)).abs() < 0.5,
            "got {}",
            fundamental.phase_degrees
        );
    }

    /// Ground truth computed by a *different* numerical method (direct Riemann-sum quadrature
    /// over the triangle wave's own closed-form definition, coded independently of this
    /// module's analytic segment formula) for an exactly-PWL unit-amplitude triangle wave over
    /// one period. This cross-checks the analytic per-segment integral against brute-force
    /// integration, not against the crate's own internal consistency.
    fn triangle(t: f64) -> f64 {
        // Period 1: 0 -> 1 over [0, 0.25], 1 -> -1 over [0.25, 0.75], -1 -> 0 over [0.75, 1].
        if t <= 0.25 {
            4.0 * t
        } else if t <= 0.75 {
            1.0 - 4.0 * (t - 0.25)
        } else {
            -1.0 + 4.0 * (t - 0.75)
        }
    }

    fn quadrature_coefficient(n: usize, steps: usize) -> Complex {
        let omega = 2.0 * PI * n as f64;
        let dt = 1.0 / steps as f64;
        let mut acc = Complex::ZERO;
        for i in 0..steps {
            let t = (i as f64 + 0.5) * dt; // midpoint rule
            let v = triangle(t);
            let (s, c) = (omega * t).sin_cos();
            acc = acc.add(
                Complex {
                    re: v * c,
                    im: -v * s,
                }
                .scale(dt),
            );
        }
        acc
    }

    #[test]
    fn four_matches_independent_quadrature_on_triangle_wave() {
        let samples = [
            Sample::new(0.0, 0.0),
            Sample::new(0.25, 1.0),
            Sample::new(0.75, -1.0),
            Sample::new(1.0, 0.0),
        ];
        let result = four(&samples, &Window::FULL, 1.0, 4).unwrap();

        for n in 1..=3usize {
            let ground_truth = quadrature_coefficient(n, 2_000_000);
            let expected_mag = 2.0 * ground_truth.abs();
            let got_mag = result.harmonics[n].magnitude;
            assert!(
                (got_mag - expected_mag).abs() < 1e-4,
                "n={n}: got {got_mag}, expected {expected_mag}"
            );
            if expected_mag > 1e-6 {
                let expected_phase = ground_truth.im.atan2(ground_truth.re).to_degrees();
                let got_phase = result.harmonics[n].phase_degrees;
                assert!(
                    (got_phase - expected_phase).abs() < 0.1,
                    "n={n}: got phase {got_phase}, expected {expected_phase}"
                );
            }
        }
        // A symmetric triangle wave has no even harmonics: n=2 magnitude must be ~0, by hand.
        assert!(
            result.harmonics[2].magnitude < 1e-9,
            "got {}",
            result.harmonics[2].magnitude
        );
    }

    /// A signal that is exactly its own fundamental sinusoid plus a known second harmonic of
    /// known relative amplitude: v(t) = sin(wt) + 0.1*sin(2wt). THD should be ~10%, by hand
    /// (ratio of the added harmonic's amplitude to the fundamental's).
    #[test]
    fn thd_of_known_two_tone_signal() {
        let freq = 10.0;
        let period = 1.0 / freq;
        let n_points = 4000;
        let samples: Vec<Sample> = (0..=n_points)
            .map(|i| {
                let t = period * (i as f64) / (n_points as f64);
                let w = 2.0 * PI * freq * t;
                Sample::new(t, w.sin() + 0.1 * (2.0 * w).sin())
            })
            .collect();
        let result = four(&samples, &Window::FULL, freq, 4).unwrap();
        assert!(
            (result.thd_percent - 10.0).abs() < 0.2,
            "got {}",
            result.thd_percent
        );
    }
}
