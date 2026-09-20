# Changelog

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
