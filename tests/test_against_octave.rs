use nalgebra;

#[cfg(test)]
mod full_solver_tests {

    #[test]
    fn silly_simple() {
        let settings = topopt::Settings::new(2, 2, 0.5)
            .with_left_bc(true, false)
            .with_bottom_right_bc(false, true)
            .with_top_left_load(0.0, -1.0)
            .with_penalty_weight(3.0)
            .with_filter_radius(1.5);
        let x = topopt::solve(settings);
        let y = topopt::mocktave::mocktave_top(2, 2, 0.5, 3.0, 1.5);
        assert!(x.relative_eq(&y, 1e-10, 0.0));
    }

    #[test]
    fn baby_beam() {
        let settings = topopt::Settings::new(10, 10, 0.5)
            .with_left_bc(true, false)
            .with_bottom_right_bc(false, true)
            .with_top_left_load(0.0, -1.0)
            .with_penalty_weight(3.0)
            .with_filter_radius(1.5);
        let x = topopt::solve(settings);
        let y = topopt::mocktave::mocktave_top(10, 10, 0.5, 3.0, 1.5);
        assert!(x.relative_eq(&y, 1e-10, 0.0));
    }

    #[test]
    fn mbb_half_resolution() {
        let settings = topopt::Settings::new(30, 10, 0.5)
            .with_left_bc(true, false)
            .with_bottom_right_bc(false, true)
            .with_top_left_load(0.0, -1.0)
            .with_penalty_weight(3.0)
            .with_filter_radius(1.5);
        let x = topopt::solve(settings);
        let y = topopt::mocktave::mocktave_top(30, 10, 0.5, 3.0, 1.5);
        assert!(x.relative_eq(&y, 1e-10, 0.0));
    }
}
