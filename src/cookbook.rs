//! # The Cookbook: Start Here for Examples, Demos, and Tips
//! # Beaucoup Boundary Conditions and Lotsa Loads
//! ## Messerschmitt–Bölkow–Blohm Simply Supported Beam (with symmetry)
//! The MBB beam is simple demo case used for topology optimization and consists of a
//! beam supported by pin joints at either end with a load in the middle. This setup
//! solves it with a symmetry constraint.
//! ```rust,no_run
#![doc = include_str!("../examples/mbb.rs")]
//! ```
//! ##  Messerschmitt–Bölkow–Blohm Simply Supported Beam (without symmetry)
//! The MBB beam is simple demo case used for topology optimization and consists of a
//! beam supported by pin joints at either end with a load in the middle. This setup
//! solves it without a symmetry constraint.
//! ```rust,no_run
#![doc = include_str!("../examples/mbb_without_reflection.rs")]
//! ```
//! # Active Elements
//! This library also supports the addition of active elements, or elements that must always
//! have material. This example shows how to define the active area.
//! ```rust,no_run
#![doc = include_str!("../examples/mbb_active.rs")]
//! ```
//! # Passive Elements
//! This library also supports the addition of active elements, or elements that must never
//! have material. This example shows how to define the passive area.
//! ```rust,no_run
#![doc = include_str!("../examples/mbb_passive.rs")]
//! ```
//!
