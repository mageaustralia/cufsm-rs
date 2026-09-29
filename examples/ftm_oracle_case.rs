//! The FTM load factors for the CalculiX oracle cases (`oracle/ftm_calculix`).
//! Prints one `case lambda` line per case, for `run_oracle.py` to compare against ccx.
use cufsm::ftm::{default_long_terms, ftm_buckle, End, Stresses, Tube};

/// Keep every case in step with `oracle/ftm_calculix/run_oracle.py`.
fn cases() -> Vec<(&'static str, Tube, Stresses)> {
    let e = 203e3;
    let nu = 0.3;
    let tube = |r: f64, t: f64, l: f64| Tube { r, t, l, e, nu };
    vec![
        (
            "long_N",
            tube(100.0, 2.0, 8000.0),
            Stresses {
                n: 1.0,
                ..Stresses::default()
            },
        ),
        (
            "med_N",
            tube(100.0, 2.0, 400.0),
            Stresses {
                n: 1.0,
                ..Stresses::default()
            },
        ),
        (
            "med_M",
            tube(100.0, 2.0, 400.0),
            Stresses {
                m: 1.0,
                ..Stresses::default()
            },
        ),
        (
            "thin_N",
            tube(500.0, 1.0, 2000.0),
            Stresses {
                n: 1.0,
                ..Stresses::default()
            },
        ),
    ]
}

fn main() {
    for (name, tube, stresses) in cases() {
        let js = default_long_terms(&tube, 6);
        let res = ftm_buckle(&tube, &stresses, End::Pinned, End::Pinned, 16, &js, 1)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        println!("{name} {:.6}", res.modes[0].load_factor);
    }
}
