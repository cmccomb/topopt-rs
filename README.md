[![build](https://github.com/cmccomb/topopt-rs/actions/workflows/tests.yml/badge.svg)](https://github.com/cmccomb/topopt-rs/actions/workflows/tests.yml)

# About
This package performs 2-dimensional topology optimization and started as a clone of ["A 99 line topology optimization code written in Matlab"](https://www.topopt.mek.dtu.dk/apps-and-software/a-99-line-topology-optimization-code-written-in-matlab)

# Usage
Usage follows exactly the same format as the original topology optimization code:
```rust
let nelx = 30;
let nely = 10;
let volfrac = 0.5;
let penalty = 3.0;
let rmin = 1.5;
let x = topopt::top(nelx, nely, volfrac, penalty, rmin);
```
