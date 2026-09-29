//! A bare shell edge, the free tip of a cantilever tube, buckles locally at
//! about half the classical stress (`ftm::End::Free`). A free end is only sound
//! opposite a clamped one: every other pairing with it is a mechanism.

use cufsm::ftm::{default_long_terms, ftm_buckle, End, Stresses, Tube};

#[test]
fn a_free_edge_buckles_at_about_half_the_classical_stress() {
    let tube = Tube {
        r: 500.0,
        t: 1.0,
        l: 2000.0,
        e: 203e3,
        nu: 0.3,
    }; // R/t = 500
    let classical = 0.605 * tube.e * tube.t / tube.r;
    let js = default_long_terms(&tube, 6);
    let n = Stresses {
        n: 1.0,
        ..Stresses::default()
    };
    let cantilever = ftm_buckle(&tube, &n, End::Clamped, End::Free, 16, &js, 1).unwrap();
    let ratio = cantilever.modes[0].load_factor / classical;
    assert!(
        (0.35..0.65).contains(&ratio),
        "free/classical = {ratio:.3}, expected about half"
    );
    // the same tube with supported ends sits on the classical stress
    let pinned = ftm_buckle(&tube, &n, End::Pinned, End::Pinned, 16, &js, 1).unwrap();
    assert!((pinned.modes[0].load_factor / classical - 1.0).abs() < 0.02);
}

#[test]
fn a_free_end_needs_a_clamped_one() {
    let tube = Tube {
        r: 500.0,
        t: 1.0,
        l: 2000.0,
        e: 203e3,
        nu: 0.3,
    };
    let js = default_long_terms(&tube, 6);
    let n = Stresses {
        n: 1.0,
        ..Stresses::default()
    };
    for (base, top) in [
        (End::Free, End::Free),
        (End::Pinned, End::Free),
        (End::Ring, End::Free),
    ] {
        let e = ftm_buckle(&tube, &n, base, top, 16, &js, 1).unwrap_err();
        assert!(
            format!("{e:?}").contains("mechanism"),
            "{base:?}-{top:?}: {e:?}"
        );
    }
}
