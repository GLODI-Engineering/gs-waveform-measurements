//! The three examples in `README.md`'s API tour, compiled and checked against values derived by
//! hand rather than recorded from a run.
//!
//! `README.md` is not included as crate documentation, so its code blocks are not compiled by
//! `cargo test --doc`. This file is the parity test that keeps them honest: if an example here
//! has to change, the README's copy of it changes in the same commit.

use gs_waveform_measurements::{
    error_norm, four, trig_targ, Edge, EventSpec, Norm, Occurrence, Sample, Window,
};

/// 10 %-to-90 % rise time of a step, both levels taken from the signal's own peak.
///
/// By hand, on the two straddling segments: the 1 us -> 1.4 us segment rises 0 -> 0.8, a slope of
/// 2e6 V/s, so it reaches 0.1 at 1.05 us. The 1.4 us -> 2 us segment rises 0.8 -> 1.0, a slope of
/// 1/3e6 V/s, so it reaches 0.9 at 1.7 us. The rise time is 650 ns, and neither level falls on a
/// sample.
#[test]
fn readme_rise_time_between_two_fractions_of_the_peak() {
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

    let t_rise = trig_targ(&at_fraction(0.1), &at_fraction(0.9)).unwrap();

    assert!(
        (t_rise - 650e-9).abs() < 1e-15,
        "rise time {t_rise} should be 650 ns"
    );
}

/// `FOUR` over one period of a symmetric triangle wave, which is exactly piecewise-linear and so
/// is integrated exactly by this crate's closed-form segment integral.
///
/// A unit triangle's harmonics are $A_n = 8/(\pi^2 n^2)$ for odd $n$ and zero for even $n$. With
/// five harmonics requested (DC plus orders 1..4) the only overtone inside the bound is the
/// third, so THD is $100\% \times A_3/A_1 = 100/9 = 11.11\%$ exactly.
#[test]
fn readme_four_of_a_triangle_matches_its_analytic_harmonics() {
    let period = 0.02; // 50 Hz
    let triangle = [
        Sample::new(0.0, 0.0),
        Sample::new(period / 4.0, 1.0),
        Sample::new(3.0 * period / 4.0, -1.0),
        Sample::new(period, 0.0),
    ];

    let spectrum = four(&triangle, &Window::new(0.0, period), 50.0, 5).unwrap();

    let a1 = spectrum.harmonics[1].magnitude;
    let a3 = spectrum.harmonics[3].magnitude;
    let analytic = |n: f64| 8.0 / (std::f64::consts::PI.powi(2) * n * n);

    assert!(
        (a1 - analytic(1.0)).abs() < 1e-12,
        "fundamental {a1} should be 8/pi^2 = {}",
        analytic(1.0)
    );
    assert!(
        (a3 - analytic(3.0)).abs() < 1e-12,
        "third harmonic {a3} should be 8/(9 pi^2) = {}",
        analytic(3.0)
    );
    assert!(
        (spectrum.harmonics[2].magnitude).abs() < 1e-12,
        "a symmetric triangle has no even harmonics"
    );
    assert!(
        (spectrum.thd_percent - 100.0 / 9.0).abs() < 1e-9,
        "THD {} should be 100/9 %",
        spectrum.thd_percent
    );
}

/// `ERROR`: the L2 norm between a measured series and a reference, evaluated at the measured
/// series' own sample times.
///
/// By hand: the two series differ only at t = 1, by 0.5, so the RMS difference over the three
/// sample times is $\sqrt{0.25/3} = 0.2887$.
#[test]
fn readme_error_norm_against_a_reference_waveform() {
    let measured = [
        Sample::new(0.0, 0.0),
        Sample::new(1.0, 1.0),
        Sample::new(2.0, 2.0),
    ];
    let reference = [
        Sample::new(0.0, 0.0),
        Sample::new(1.0, 1.5),
        Sample::new(2.0, 2.0),
    ];

    let l2 = error_norm(&measured, &reference, Norm::L2).unwrap();

    assert!(
        (l2 - (0.25f64 / 3.0).sqrt()).abs() < 1e-15,
        "L2 norm {l2} should be sqrt(0.25/3)"
    );
}
