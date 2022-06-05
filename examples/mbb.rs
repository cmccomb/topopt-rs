fn main() {
    let settings = topopt::Settings::new(60, 20, 0.5)
        .with_left_bc(true, false)
        .with_bottom_right_bc(false, true)
        .with_top_left_load(0.0, -1.0);

    topopt::solve(settings);
}
