[![build](https://github.com/cmccomb/topopt-rs/actions/workflows/tests.yml/badge.svg)](https://github.com/cmccomb/topopt-rs/actions/workflows/tests.yml)
[![Crates.io](https://img.shields.io/crates/v/topopt.svg)](https://crates.io/crates/topopt)
[![docs.rs](https://docs.rs/topopt/badge.svg)](https://docs.rs/topopt)

# About
This package performs 2-dimensional topology optimization and is a port of ["A 99 line topology optimization code written in Matlab"](https://www.topopt.mek.dtu.dk/apps-and-software/a-99-line-topology-optimization-code-written-in-matlab).

# Usage
Usage follows almost exactly the same format as the original topology optimization code:
```rust
let nelx = 60;
let nely = 10;
let volfrac = 0.5;
let penalty = 3.0;
let rmin = 1.5;
let x = topopt::top(nelx, nely, volfrac, penalty, rmin, None, None, None, None);
```
This will display the progress of the algorithm and a visualization of the optimized structure in the command line

![](https://raw.githubusercontent.com/cmccomb/topopt-rs/master/mbb.gif)