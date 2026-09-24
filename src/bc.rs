//! The five longitudinal integrals, CUFSM `BC_I1_5.m`.
//!
//! With `Ym` the longitudinal shape function of term `m` for the chosen boundary condition, over
//! `0..a`:
//!
//! - `I1 = ∫ Ym Yn`
//! - `I2 = ∫ Ym'' Yn`
//! - `I3 = ∫ Ym Yn''`
//! - `I4 = ∫ Ym'' Yn''`
//! - `I5 = ∫ Ym' Yn'`
//!
//! Ported expression for expression, so a term-by-term comparison with CUFSM holds to round-off.

use crate::model::BoundaryCondition;
use std::f64::consts::PI;

/// `(-1)^(k-1)` for a whole number `k`.
fn alt(k: f64) -> f64 {
    if (k - 1.0).rem_euclid(2.0) == 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// `[I1, I2, I3, I4, I5]` for longitudinal terms `kk` and `nn` over length `a`.
pub fn bc_i1_5(bc: BoundaryCondition, kk: f64, nn: f64, a: f64) -> [f64; 5] {
    let pi2 = PI * PI;
    let pi4 = pi2 * pi2;
    match bc {
        BoundaryCondition::SS => {
            if kk == nn {
                [
                    a / 2.0,
                    -kk * kk * pi2 / a / 2.0,
                    -nn * nn * pi2 / a / 2.0,
                    pi4 * kk.powi(4) / 2.0 / a.powi(3),
                    pi2 * kk * kk / 2.0 / a,
                ]
            } else {
                [0.0; 5]
            }
        }
        BoundaryCondition::CC => {
            if kk == nn {
                let i1 = if kk == 1.0 { 3.0 * a / 8.0 } else { a / 4.0 };
                [
                    i1,
                    -(kk * kk + 1.0) * pi2 / 4.0 / a,
                    -(nn * nn + 1.0) * pi2 / 4.0 / a,
                    pi4 * ((kk * kk + 1.0).powi(2) + 4.0 * kk * kk) / 4.0 / a.powi(3),
                    (1.0 + kk * kk) * pi2 / 4.0 / a,
                ]
            } else if kk - nn == 2.0 {
                [
                    -a / 8.0,
                    (kk * kk + 1.0) * pi2 / 8.0 / a - kk * pi2 / 4.0 / a,
                    (nn * nn + 1.0) * pi2 / 8.0 / a + nn * pi2 / 4.0 / a,
                    -(kk - 1.0).powi(2) * (nn + 1.0).powi(2) * pi4 / 8.0 / a.powi(3),
                    -(1.0 + kk * nn) * pi2 / 8.0 / a,
                ]
            } else if kk - nn == -2.0 {
                [
                    -a / 8.0,
                    (kk * kk + 1.0) * pi2 / 8.0 / a + kk * pi2 / 4.0 / a,
                    (nn * nn + 1.0) * pi2 / 8.0 / a - nn * pi2 / 4.0 / a,
                    -(kk + 1.0).powi(2) * (nn - 1.0).powi(2) * pi4 / 8.0 / a.powi(3),
                    -(1.0 + kk * nn) * pi2 / 8.0 / a,
                ]
            } else {
                [0.0; 5]
            }
        }
        BoundaryCondition::SC => {
            if kk == nn {
                [
                    (1.0 + (kk + 1.0).powi(2) / (kk * kk)) * a / 2.0,
                    -(kk + 1.0).powi(2) * pi2 / a,
                    -(kk + 1.0).powi(2) * pi2 / a,
                    (kk + 1.0).powi(2) * pi4 * ((kk + 1.0).powi(2) + kk * kk) / 2.0 / a.powi(3),
                    (1.0 + kk).powi(2) * pi2 / a,
                ]
            } else if kk - nn == 1.0 {
                [
                    (kk + 1.0) * a / 2.0 / kk,
                    -(kk + 1.0) * kk * pi2 / 2.0 / a,
                    -(nn + 1.0).powi(2) * pi2 * (kk + 1.0) / 2.0 / a / kk,
                    (kk + 1.0) * kk * (nn + 1.0).powi(2) * pi4 / 2.0 / a.powi(3),
                    (1.0 + kk) * (1.0 + nn) * pi2 / 2.0 / a,
                ]
            } else if kk - nn == -1.0 {
                [
                    (nn + 1.0) * a / 2.0 / nn,
                    -(kk + 1.0).powi(2) * pi2 * (nn + 1.0) / 2.0 / a / nn,
                    -(nn + 1.0) * nn * pi2 / 2.0 / a,
                    (kk + 1.0).powi(2) * nn * (nn + 1.0) * pi4 / 2.0 / a.powi(3),
                    (1.0 + kk) * (1.0 + nn) * pi2 / 2.0 / a,
                ]
            } else {
                [0.0; 5]
            }
        }
        BoundaryCondition::CF => {
            let hk = kk - 0.5;
            let hn = nn - 0.5;
            if kk == nn {
                [
                    3.0 * a / 2.0 - 2.0 * a * alt(kk) / hk / PI,
                    hk * hk * pi2 * (alt(kk) / hk / PI - 0.5) / a,
                    hn * hn * pi2 * (alt(nn) / hn / PI - 0.5) / a,
                    hk.powi(4) * pi4 / 2.0 / a.powi(3),
                    hk * hk * pi2 / 2.0 / a,
                ]
            } else {
                [
                    a - a * alt(kk) / hk / PI - a * alt(nn) / hn / PI,
                    hk * hk * pi2 * (alt(kk) / hk / PI) / a,
                    hn * hn * pi2 * (alt(nn) / hn / PI) / a,
                    0.0,
                    0.0,
                ]
            }
        }
        BoundaryCondition::CG => {
            let hk = kk - 0.5;
            let hn = nn - 0.5;
            if kk == nn {
                let i1 = if kk == 1.0 { 3.0 * a / 8.0 } else { a / 4.0 };
                [
                    i1,
                    -(hk * hk + 0.25) * pi2 / a / 4.0,
                    -(hk * hk + 0.25) * pi2 / a / 4.0,
                    (hk * hk + 0.25).powi(2) * pi4 / 4.0 / a.powi(3)
                        + hk * hk * pi4 / 4.0 / a.powi(3),
                    hk * hk * pi2 / a / 4.0 + pi2 / 16.0 / a,
                ]
            } else if kk - nn == 1.0 {
                [
                    -a / 8.0,
                    (hk * hk + 0.25) * pi2 / a / 8.0 - hk * pi2 / a / 8.0,
                    (hn * hn + 0.25) * pi2 / a / 8.0 + hn * pi2 / a / 8.0,
                    -nn.powi(4) * pi4 / 8.0 / a.powi(3),
                    -nn * nn * pi2 / 8.0 / a,
                ]
            } else if kk - nn == -1.0 {
                [
                    -a / 8.0,
                    (hk * hk + 0.25) * pi2 / a / 8.0 + hk * pi2 / a / 8.0,
                    (hn * hn + 0.25) * pi2 / a / 8.0 - hn * pi2 / a / 8.0,
                    -kk.powi(4) * pi4 / 8.0 / a.powi(3),
                    -kk * kk * pi2 / 8.0 / a,
                ]
            } else {
                [0.0; 5]
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// I1..I5 for S-S are the textbook integrals of sin(m pi y / a).
    #[test]
    fn simply_supported_integrals_are_the_sine_integrals() {
        let a = 3.0;
        let [i1, i2, i3, i4, i5] = bc_i1_5(BoundaryCondition::SS, 2.0, 2.0, a);
        let c = 2.0 * PI / a;
        assert!((i1 - a / 2.0).abs() < 1e-15);
        assert!((i2 + c * c * a / 2.0).abs() < 1e-12);
        assert!((i3 - i2).abs() < 1e-15);
        assert!((i4 - c.powi(4) * a / 2.0).abs() < 1e-10);
        assert!((i5 - c * c * a / 2.0).abs() < 1e-12);
        assert_eq!(bc_i1_5(BoundaryCondition::SS, 1.0, 2.0, a), [0.0; 5]);
    }

    #[test]
    fn alternating_sign() {
        assert_eq!(alt(1.0), 1.0);
        assert_eq!(alt(2.0), -1.0);
        assert_eq!(alt(3.0), 1.0);
    }
}
