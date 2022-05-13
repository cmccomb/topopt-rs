fn main() {
    // Define a mask of
    let mask = nalgebra::DMatrix::from_fn(20, 60, |idx, jdx| {
        if ((idx as f64 - 9.5).powf(2.0) + (jdx as f64 - 29.5).powf(2.0)).sqrt() < 8.0 {
            true
        } else {
            false
        }
    });
    topopt::solve(topopt::Settings::default().with_active_elements(mask));
}
