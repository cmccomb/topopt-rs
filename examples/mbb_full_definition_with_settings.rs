use nalgebra::DMatrix;
use topopt::{solve, Settings};

fn main() {
    let mut loads = DMatrix::from_element(121, 31, (0.0, 0.0));
    loads[(60, 0)] = (0.0, -1.0);

    let settings = Settings::new(120, 30, 0.5)
        .with_bottom_right_boundary(false, true)
        .with_bottom_left_boundary(false, true)
        .with_vertical_midline_boundary(true, false)
        .with_loads(loads);

    // Solve
    solve(settings);
}
