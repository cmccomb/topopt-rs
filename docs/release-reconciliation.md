# Release reconciliation

The crates.io `topopt 0.1.7` archive was published from a dirty working tree based
on commit `74ba553`. Its VCS metadata identifies that base, not the full source.
The recovered branch `codex/published-0.1.7`, commit
`530ae3d`, records the actual published source without changing that historical API.

Archive: https://static.crates.io/crates/topopt/topopt-0.1.7.crate

SHA-256: `56e4fee7956081af9c2399043a34ea4bff0939118537105bdcf228e57ab36df8`

All 16 authored files included in the archive match the recovered branch byte for
byte, treating `Cargo.toml.orig` as the source manifest. Cargo-generated
`Cargo.toml`, `Cargo.lock`, and `.cargo_vcs_info.json` are excluded from that source
comparison. No release tag was rewritten.

## Integrated 0.2.0 candidate

The candidate combines the released mask validation/output behavior with current
mainline finite-element load assembly and right-boundary fixes. It removes the
local mocktave stub, selects mocktave 0.1.6's native backend, and bounds the matrix
libraries to a compatible pair. Native Octave is needed only for the optional
reference tests: `cargo test --all-features --release`.

`set_load` now accepts numerical force components. Existing code passing booleans
must migrate: call a boundary-condition setter when constraining nodes, or pass
floating-point forces when applying loads. This API correction requires 0.2.0.

Validation compares element stiffness and four optimized density fields with the
original Octave formulation, from a 2-by-2 grid through a 60-by-20 MBB beam. Mask,
load, boundary, and terminal-output regression tests also run. These are numerical
regression checks, not a general certification of every boundary configuration.
