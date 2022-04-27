fn main() {
    let x = topopt::top(30, 10, 0.5, 3.0, 1.5);
    let (nx, ny) = x.shape();
    for ex in 0..nx {
        for ey in 0..ny {
            if x[(ex, ey)] > 0.5 {
                print!("█");
            } else {
                print!(" ");
            }
        }
        print!("\n");
    }
}
