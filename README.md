[![build](https://github.com/cmccomb/topopt-rs/actions/workflows/tests.yml/badge.svg)](https://github.com/cmccomb/topopt-rs/actions/workflows/tests.yml)
[![Crates.io](https://img.shields.io/crates/v/topopt.svg)](https://crates.io/crates/topopt)
[![docs.rs](https://docs.rs/topopt/badge.svg)](https://docs.rs/topopt)

# About
This package performs 2-dimensional topology optimization and is a port of ["A 99 line topology optimization code written in Matlab"](https://www.topopt.mek.dtu.dk/apps-and-software/a-99-line-topology-optimization-code-written-in-matlab).

# Basic Usage
Running the solve function with default settings will find a solution to the Messerschmitt–Bölkow–Blohm simply supported beam (enforcing symmetry).
```rust
topopt::solve(topopt::Settings::default());
```
In an interactive terminal, each iteration redraws the density visualization and
shows aligned `Iteration`, `Objective`, `Volume`, and `Change` values. Redrawing
clears the display without resetting terminal settings. When stdout is redirected
to a file or pipe, or `TERM=dumb`, output is one plain-text progress line per
iteration, without a visualization or terminal control codes:

```text
Iteration:    1  Objective:    123.457  Volume: 0.500  Change: 0.123
```

![](https://raw.githubusercontent.com/cmccomb/topopt-rs/master/mbb.gif)

Alternatively, we could set up with the same simulation explicitly:
```rust
topopt::solve(
    topopt::Settings::new(60, 20, 0.5)
        .with_left_bc(true, false)
        .with_bottom_right_bc(false, true)
        .with_top_left_load(0.0, -1.0),
);
```
