use nalgebra::DMatrix;

fn main() {
    // // Define MBB load case
    let mut loads = DMatrix::from_element(61, 21, (0.0, 0.0));
    loads[(0, 0)] = (0.0, -1.0);

    // // Define MBB boundary condition
    let boundary = DMatrix::from_fn(61, 21, |idx, jdx| {
        if idx == 0 {
            (true, false)
        } else if idx == 60 && jdx == 20 {
            (false, true)
        } else {
            (false, false)
        }
    });

    // Solve
    topopt::top(
        60,
        20,
        0.5,
        3.0,
        1.5,
        Some(loads),
        Some(boundary),
        None,
        None,
    );
}
