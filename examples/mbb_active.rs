use nalgebra::DMatrix;

fn main() {
    let m = DMatrix::from_fn(10, 30, |idx, jdx| {
        if ((idx as f64 - 4.5).powf(2.0) + (jdx as f64 - 14.5).powf(2.0)).sqrt() < 4.0 {
            true
        } else {
            false
        }
    });
    println!("{m}");
    topopt::top(30, 10, 0.5, 3.0, 1.5, None, None, None, Some(m));
}
