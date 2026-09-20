# Changelog

## 0.2.0

- Reconcile the exact published 0.1.7 mask validation and terminal output with
  subsequent mainline load assembly and right-boundary fixes.
- **Breaking:** `set_load(idx, jdx, x, y)` now takes floating-point force
  components. In 0.1.7 the boolean arguments incorrectly changed supports;
  use the boundary-condition setters to constrain nodes.
- Bound nalgebra/nalgebra-sparse to compatible 0.35/0.12 minor versions.
- Remove the placeholder mocktave patch and run the real native Octave reference
  comparisons, including the full 60-by-20 MBB beam, in CI.
- Record published archive provenance in `docs/release-reconciliation.md`.

## 0.1.7 (2026-09-19)

- Validate incoming active and passive element masks against `(nely, nelx)`
  before changing settings. Invalid masks now panic at the setter with the
  expected and actual dimensions.
- Correct the passive-element documentation to describe void elements.
- Add rectangular-grid regression coverage for mask dimensions, element
  positions, and preserving settings when a mask is rejected.
- Replace the full terminal reset with a screen redraw and aligned progress
  labels. Redirected stdout and `TERM=dumb` use one plain-text progress line per
  iteration, without animation or terminal control codes.

This patch is based on the published 0.1.6 release and preserves its public
method signatures. Later changes on `master` are outside this release.
