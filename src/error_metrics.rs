//! Waveform-comparison measurements: `ERR`/`ERR1`/`ERR2` and `ERROR`.
//!
//! These are two genuinely different computations that Xyce's own reference guide is careful to
//! distinguish ("the ERROR measure is not the functional equivalent of the ERR1 or ERR2
//! measures"): `ERR1`/`ERR2` are a per-point relative-difference statistic between two
//! *simulation* variables (typically compared at the same time points), while `ERROR` is a norm
//! (`L1`/`L2`/`INFNORM`) between a measured waveform and a separate reference waveform,
//! interpolated onto a common time base. Both fit this crate's "two time series in, a number
//! out" shape; reading Xyce's `FILE=` reference-waveform format itself is left to a later
//! integration layer (see this crate's own top-level docs).

use crate::sample::{interpolate, MeasureError, Result, Sample};

/// `ERR1`: the RMS of the per-sample relative difference between `measured` and `comparison`.
///
/// $$\mathrm{ERR1} = \sqrt{\frac{1}{N}\sum_{i=1}^N \left(\frac{M_i -
/// C_i}{\max(\mathrm{minval}, |M_i|)}\right)^2}$$
///
/// `comparison` is linearly interpolated onto `measured`'s own sample times (so the two series
/// need not share a time grid); a sample is skipped if `|M_i|` falls outside `[ymin, ymax]`.
pub fn err1(
    measured: &[Sample],
    comparison: &[Sample],
    minval: f64,
    ymin: f64,
    ymax: f64,
) -> Result<f64> {
    let terms = err_terms(measured, comparison, minval, ymin, ymax)?;
    let n = terms.len() as f64;
    let sum_sq: f64 = terms.iter().map(|t| t * t).sum();
    Ok((sum_sq / n).sqrt())
}

/// `ERR2`: the mean of the per-sample relative difference between `measured` and `comparison`
/// (same terms as [`err1`], but averaged directly instead of RMS'd).
///
/// $$\mathrm{ERR2} = \frac{1}{N}\sum_{i=1}^N \frac{M_i - C_i}{\max(\mathrm{minval}, |M_i|)}$$
pub fn err2(
    measured: &[Sample],
    comparison: &[Sample],
    minval: f64,
    ymin: f64,
    ymax: f64,
) -> Result<f64> {
    let terms = err_terms(measured, comparison, minval, ymin, ymax)?;
    let n = terms.len() as f64;
    Ok(terms.iter().sum::<f64>() / n)
}

fn err_terms(
    measured: &[Sample],
    comparison: &[Sample],
    minval: f64,
    ymin: f64,
    ymax: f64,
) -> Result<Vec<f64>> {
    if measured.is_empty() || comparison.is_empty() {
        return Err(MeasureError::EmptySeries);
    }
    let terms: Vec<f64> = measured
        .iter()
        .filter(|m| m.v.abs() >= ymin && m.v.abs() <= ymax)
        .filter_map(|m| {
            let c = interpolate(comparison, m.t)?;
            let denom = minval.max(m.v.abs());
            Some((m.v - c) / denom)
        })
        .collect();
    if terms.is_empty() {
        return Err(MeasureError::EmptyWindow);
    }
    Ok(terms)
}

/// Which norm [`error_norm`] computes — Xyce's `COMP_FUNCTION=`/`ERROR` measure norms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Norm {
    /// Mean absolute difference: $\frac{1}{N}\sum_i |M_i - R_i|$.
    L1,
    /// Root-mean-square difference: $\sqrt{\frac{1}{N}\sum_i (M_i - R_i)^2}$ (the default, per
    /// Xyce's own `ERROR` measure).
    L2,
    /// Maximum absolute difference: $\max_i |M_i - R_i|$.
    InfNorm,
}

/// `ERROR`: the norm between `measured` and a `reference` waveform, evaluated at `measured`'s
/// own sample times restricted to the overlap of both series' domains (`reference` is linearly
/// interpolated onto that grid) — the norm computation itself, over two given time series.
/// Reading a reference waveform out of a file (Xyce's `FILE=`/`INDEPVARCOL=`/`DEPVARCOL=`
/// qualifiers) is file-I/O and format-parsing left to a later integration layer.
pub fn error_norm(measured: &[Sample], reference: &[Sample], norm: Norm) -> Result<f64> {
    if measured.is_empty() || reference.is_empty() {
        return Err(MeasureError::EmptySeries);
    }
    let lo = reference[0].t;
    let hi = reference[reference.len() - 1].t;
    let diffs: Vec<f64> = measured
        .iter()
        .filter(|m| m.t >= lo && m.t <= hi)
        .filter_map(|m| interpolate(reference, m.t).map(|r| m.v - r))
        .collect();
    if diffs.is_empty() {
        return Err(MeasureError::EmptyWindow);
    }
    let n = diffs.len() as f64;
    Ok(match norm {
        Norm::L1 => diffs.iter().map(|d| d.abs()).sum::<f64>() / n,
        Norm::L2 => (diffs.iter().map(|d| d * d).sum::<f64>() / n).sqrt(),
        Norm::InfNorm => diffs.iter().fold(0.0_f64, |acc, d| acc.max(d.abs())),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// measured = [1,2,3], comparison identical: ERR1 and ERR2 must be exactly 0, by hand.
    #[test]
    fn err_of_identical_series_is_zero() {
        let m = [
            Sample::new(0.0, 1.0),
            Sample::new(1.0, 2.0),
            Sample::new(2.0, 3.0),
        ];
        assert!(err1(&m, &m, 1e-12, 1e-15, 1e15).unwrap().abs() < 1e-12);
        assert!(err2(&m, &m, 1e-12, 1e-15, 1e15).unwrap().abs() < 1e-12);
    }

    /// measured = [2,4], comparison = [1,2] (comparison is always half): each term = (M-C)/M =
    /// 0.5 exactly. ERR2 = mean(0.5,0.5) = 0.5. ERR1 = sqrt(mean(0.25,0.25)) = 0.5. By hand.
    #[test]
    fn err1_err2_hand_computed() {
        let m = [Sample::new(0.0, 2.0), Sample::new(1.0, 4.0)];
        let c = [Sample::new(0.0, 1.0), Sample::new(1.0, 2.0)];
        assert!((err2(&m, &c, 1e-12, 1e-15, 1e15).unwrap() - 0.5).abs() < 1e-12);
        assert!((err1(&m, &c, 1e-12, 1e-15, 1e15).unwrap() - 0.5).abs() < 1e-12);
    }

    /// measured = reference + constant offset of 3 everywhere: L1 = 3, L2 = 3, InfNorm = 3, by
    /// hand (a constant-offset difference makes all three norms coincide).
    #[test]
    fn error_norm_constant_offset() {
        let reference = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 10.0),
            Sample::new(2.0, 0.0),
        ];
        let measured = [
            Sample::new(0.0, 3.0),
            Sample::new(1.0, 13.0),
            Sample::new(2.0, 3.0),
        ];
        assert!((error_norm(&measured, &reference, Norm::L1).unwrap() - 3.0).abs() < 1e-9);
        assert!((error_norm(&measured, &reference, Norm::L2).unwrap() - 3.0).abs() < 1e-9);
        assert!((error_norm(&measured, &reference, Norm::InfNorm).unwrap() - 3.0).abs() < 1e-9);
    }

    /// A single outlier: measured = reference except one point off by 10 among four points ->
    /// L1 = 10/4 = 2.5, L2 = sqrt(100/4) = 5.0, InfNorm = 10, by hand.
    #[test]
    fn error_norm_single_outlier() {
        let reference = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 0.0),
            Sample::new(2.0, 0.0),
            Sample::new(3.0, 0.0),
        ];
        let measured = [
            Sample::new(0.0, 0.0),
            Sample::new(1.0, 10.0),
            Sample::new(2.0, 0.0),
            Sample::new(3.0, 0.0),
        ];
        assert!((error_norm(&measured, &reference, Norm::L1).unwrap() - 2.5).abs() < 1e-9);
        assert!((error_norm(&measured, &reference, Norm::L2).unwrap() - 5.0).abs() < 1e-9);
        assert!((error_norm(&measured, &reference, Norm::InfNorm).unwrap() - 10.0).abs() < 1e-9);
    }
}
