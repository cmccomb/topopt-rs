fn main() {
    let settings = topopt::Settings::new(120, 20, 0.5)
        .with_bottom_right_bc(false, true)
        .with_bottom_left_bc(false, true)
        .with_vertical_midline_bc(true, false)
        .with_top_middle_load(0.0, -1.0);

    topopt::solve(settings);
}
