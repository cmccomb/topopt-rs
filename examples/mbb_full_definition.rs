use nalgebra::DMatrix;

fn main() {
    // // Define MBB load case
    let mut loads = DMatrix::from_element(121, 31, (0.0, 0.0));
    loads[(60, 0)] = (0.0, -1.0);

    // // Define MBB boundary condition
    let boundary = DMatrix::from_fn(121, 31, |idx, jdx| {
        if idx == 60 {
            (true, false)
        } else if idx == 120 && jdx == 30 {
            (false, true)
        } else if idx == 0 && jdx == 30 {
            (false, true)
        } else {
            (false, false)
        }
    });

    // Solve
    topopt::top(
        120,
        30,
        0.5,
        3.0,
        1.5,
        Some(loads),
        Some(boundary),
        None,
        None,
    );
}
