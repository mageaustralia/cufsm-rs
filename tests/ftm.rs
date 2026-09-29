//! The finite tube method against closed-form and code results for circular tubes.

use cufsm::ftm::{default_long_terms, ftm_buckle, ftm_field, End, Stresses, Tube};
use std::f64::consts::PI;

const E: f64 = 210_000.0;
const NU: f64 = 0.3;

fn tube(r: f64, t: f64, l: f64) -> Tube {
    Tube {
        r,
        t,
        l,
        e: E,
        nu: NU,
    }
}

fn only(k: usize) -> Stresses {
    let mut s = Stresses::default();
    match k {
        0 => s.n = 1.0,
        1 => s.m = 1.0,
        2 => s.t = 1.0,
        _ => s.v = 1.0,
    }
    s
}

/// The classical axial buckling stress of a cylinder, E t / (R √(3(1 − ν²))).
fn classical(r: f64, t: f64) -> f64 {
    E * t / r / (3.0 * (1.0 - NU * NU)).sqrt()
}

fn near(got: f64, want: f64, tol: f64, what: &str) {
    assert!(
        (got / want - 1.0).abs() < tol,
        "{what}: {got} vs {want} ({:+.2} %)",
        (got / want - 1.0) * 100.0
    );
}

/// A long tube is an Euler column: π²E R² / (2 (kL)²) for every set of end conditions.
#[test]
fn a_long_tube_buckles_as_an_euler_column() {
    let (r, l) = (100.0, 20_000.0);
    let euler = |k: f64| PI * PI * E * r * r / (2.0 * (k * l).powi(2));
    let js: Vec<usize> = (1..=10).collect();
    for (base, top, k) in [
        (End::Pinned, End::Pinned, 1.0),
        (End::Clamped, End::Free, 2.0),
        (End::Clamped, End::Clamped, 0.5),
        (End::Clamped, End::Pinned, 0.699_155),
    ] {
        let r1 = ftm_buckle(&tube(r, 2.0, l), &only(0), base, top, 3, &js, 1).unwrap();
        near(
            r1.modes[0].load_factor,
            euler(k),
            0.003,
            &format!("{base:?}-{top:?}"),
        );
        assert_eq!(
            r1.modes[0].circ_waves, 1,
            "a column bends in the first harmonic"
        );
    }
}

/// A medium-length cylinder in compression buckles at the classical stress; a long one a few
/// per cent under it, as EN 1993-1-6 D.1.2.1's long-cylinder factor has it (0.967 here).
#[test]
fn compression_reaches_the_classical_stress() {
    let medium = tube(100.0, 1.0, 300.0);
    let r1 = ftm_buckle(
        &medium,
        &only(0),
        End::Pinned,
        End::Pinned,
        12,
        &default_long_terms(&medium, 8),
        1,
    )
    .unwrap();
    near(
        r1.modes[0].load_factor,
        classical(100.0, 1.0),
        0.02,
        "medium",
    );
    let long = tube(100.0, 1.0, 1000.0);
    let r2 = ftm_buckle(
        &long,
        &only(0),
        End::Pinned,
        End::Pinned,
        12,
        &default_long_terms(&long, 8),
        1,
    )
    .unwrap();
    let omega = 1000.0 / (100.0f64 * 1.0).sqrt();
    let cx = 1.0 + 0.2 / 6.0 * (1.0 - 2.0 * omega * 1.0 / 100.0);
    near(
        r2.modes[0].load_factor,
        cx * classical(100.0, 1.0),
        0.02,
        "long",
    );
}

/// Bending: the peak compressive stress at buckling is close to the classical axial stress.
#[test]
fn bending_reaches_the_classical_stress() {
    let t = tube(100.0, 1.0, 300.0);
    let r = ftm_buckle(
        &t,
        &only(1),
        End::Pinned,
        End::Pinned,
        16,
        &default_long_terms(&t, 4),
        1,
    )
    .unwrap();
    near(
        r.modes[0].load_factor,
        classical(100.0, 1.0),
        0.03,
        "bending",
    );
}

/// Torsion of a medium-length cylinder, EN 1993-1-6 D.1.4: τ = 0.75 E (t/r) √(1/ω), ω = l/√(rt).
#[test]
fn torsion_agrees_with_en_1993_1_6() {
    let t = tube(100.0, 1.0, 1000.0);
    let r = ftm_buckle(
        &t,
        &only(2),
        End::Pinned,
        End::Pinned,
        16,
        &default_long_terms(&t, 6),
        1,
    )
    .unwrap();
    let tau = 0.75 * E * (1.0 / 100.0) * (1.0 / (1000.0 / 10.0f64)).sqrt();
    near(r.modes[0].load_factor, tau, 0.08, "torsion");
}

/// The load factor on section actions: an Euler column under a compressive force N buckles at
/// λN = π²EI/L², I = πR³t.
#[test]
fn actions_give_the_euler_force() {
    let t = tube(100.0, 2.0, 20_000.0);
    let n = 10_000.0;
    let s = Stresses::from_actions(&t, n, 0.0, 0.0, 0.0);
    let r = ftm_buckle(&t, &s, End::Pinned, End::Pinned, 3, &[1, 2, 3, 4], 1).unwrap();
    let pe = PI * PI * E * PI * 100.0f64.powi(3) * 2.0 / 20_000.0f64.powi(2);
    near(r.modes[0].load_factor * n, pe, 0.003, "Euler force");
}

/// A load factor scales with E and inversely with the stress; bending one way or the other is the
/// same problem.
#[test]
fn invariants() {
    let t = tube(100.0, 1.0, 300.0);
    let js = default_long_terms(&t, 4);
    let base = ftm_buckle(&t, &only(1), End::Pinned, End::Pinned, 12, &js, 1)
        .unwrap()
        .modes[0]
        .load_factor;
    let stiff = ftm_buckle(
        &Tube { e: 2.0 * E, ..t },
        &only(1),
        End::Pinned,
        End::Pinned,
        12,
        &js,
        1,
    )
    .unwrap();
    near(stiff.modes[0].load_factor, 2.0 * base, 1e-9, "E");
    let s3 = Stresses {
        m: 3.0,
        ..Default::default()
    };
    near(
        ftm_buckle(&t, &s3, End::Pinned, End::Pinned, 12, &js, 1)
            .unwrap()
            .modes[0]
            .load_factor,
        base / 3.0,
        1e-9,
        "stress",
    );
    let neg = Stresses {
        m: -1.0,
        ..Default::default()
    };
    near(
        ftm_buckle(&t, &neg, End::Pinned, End::Pinned, 12, &js, 1)
            .unwrap()
            .modes[0]
            .load_factor,
        base,
        1e-6,
        "sign",
    );
}

/// A free-free tube, one held only by a pin, and a tube with no stress are refused.
#[test]
fn refusals() {
    let t = tube(100.0, 1.0, 300.0);
    for (b, top) in [(End::Free, End::Free), (End::Pinned, End::Free)] {
        let e = ftm_buckle(&t, &only(0), b, top, 4, &[1, 2], 1)
            .unwrap_err()
            .to_string();
        assert!(e.contains("mechanism"), "{b:?}-{top:?}: {e}");
    }
    assert!(ftm_buckle(
        &t,
        &Stresses::default(),
        End::Pinned,
        End::Pinned,
        4,
        &[1],
        1
    )
    .is_err());
    assert!(ftm_buckle(
        &tube(1.0, 2.0, 300.0),
        &only(0),
        End::Pinned,
        End::Pinned,
        4,
        &[1],
        1
    )
    .is_err());
}

/// The mode field: a pinned column's first mode is a half sine, zero at both ends, peaking mid
/// length in the first harmonic.
#[test]
fn the_mode_field_is_the_buckled_shape() {
    let t = tube(100.0, 2.0, 20_000.0);
    let r = ftm_buckle(&t, &only(0), End::Pinned, End::Pinned, 3, &[1, 2, 3, 4], 1).unwrap();
    let f = ftm_field(&r, &r.modes[0], 8, 5);
    let w = |i: usize, j: usize| f[i * 5 + j][2];
    let peak = (0..8).map(|i| w(i, 2).abs()).fold(0.0, f64::max);
    for i in 0..8 {
        assert!(
            w(i, 0).abs() < 1e-6 * peak && w(i, 4).abs() < 1e-6 * peak,
            "ends move: {} {}",
            w(i, 0),
            w(i, 4)
        );
    }
    // First harmonic: opposite sides move the same way in space, so w(θ) = −w(θ + π).
    for i in 0..4 {
        assert!((w(i, 2) + w(i + 4, 2)).abs() < 1e-6 * peak);
    }
}

/// A cantilever under a moment: with a bare free top the compression side buckles at that edge,
/// at about half the classical stress; a stiff ring there restores the classical value, and does
/// not stiffen the tube's sway as a column.
#[test]
fn a_ring_keeps_a_free_end_round() {
    let t = tube(251.0, 6.0, 6000.0);
    let js = default_long_terms(&t, 6);
    let m = |top| {
        ftm_buckle(&t, &only(1), End::Clamped, top, 16, &js, 1)
            .unwrap()
            .modes[0]
            .load_factor
    };
    let (bare, ring) = (m(End::Free), m(End::Ring));
    let cl = classical(251.0, 6.0);
    assert!(
        bare < 0.6 * cl && bare > 0.35 * cl,
        "bare edge {bare} vs classical {cl}"
    );
    near(ring, cl, 0.04, "ringed top");
    let n = |top| {
        ftm_buckle(
            &t,
            &only(0),
            End::Clamped,
            top,
            4,
            &(1..=10).collect::<Vec<_>>(),
            1,
        )
        .unwrap()
        .modes[0]
            .load_factor
    };
    near(
        n(End::Ring),
        n(End::Free),
        0.002,
        "sway with and without the ring",
    );
    let e = ftm_buckle(&t, &only(0), End::Ring, End::Ring, 4, &[1, 2], 1)
        .unwrap_err()
        .to_string();
    assert!(e.contains("mechanism"), "{e}");
}
