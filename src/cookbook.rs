//! # The Cookbook
//! # Beaucoup Boundary Conditions and Lotsa Loads
//! ## Messerschmitt–Bölkow–Blohm Simply Supported Beam (with symmetry)
//! ```rust
//! let settings = topopt::Settings::new(60, 20, 0.5)
//!     .with_left_bc(true, false)
//!     .with_bottom_right_bc(false, true)
//!     .with_top_left_load(0.0, -1.0);
//!
//! topopt::solve(settings);
//! ```
//! ##  Messerschmitt–Bölkow–Blohm Simply Supported Beam (without symmetry)
//! ```rust
//! let settings = topopt::Settings::new(120, 20, 0.5)
//!     .with_bottom_right_bc(false, true)
//!     .with_bottom_left_bc(true, true)
//!     .with_top_middle_load(0.0, -1.0);
//!
//! topopt::solve(settings);
//! ```
//! # Active Elements
//! # Passive Elements
//!
