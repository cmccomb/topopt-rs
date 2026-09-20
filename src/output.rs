//! Terminal visualization and plain-text solver progress.

use nalgebra::DMatrix;
use std::io::{self, IsTerminal, Write};

/// Print one iteration, animating only when stdout is an interactive terminal.
pub(crate) fn print_progress(
    iteration: usize,
    objective: f64,
    volume: f64,
    change: f64,
    densities: &DMatrix<f64>,
) {
    let interactive = io::stdout().is_terminal()
        && !matches!(std::env::var_os("TERM"), Some(term) if term == "dumb");
    let frame = format_progress(iteration, objective, volume, change, densities, interactive);
    print!("{frame}");
    io::stdout()
        .flush()
        .expect("Could not flush solver progress");
}

/// Format a complete frame so it can be printed in one operation.
fn format_progress(
    iteration: usize,
    objective: f64,
    volume: f64,
    change: f64,
    densities: &DMatrix<f64>,
    interactive: bool,
) -> String {
    let mut output = String::new();
    if interactive {
        // Clear the display and return to the top without resetting the terminal.
        output.push_str("\x1b[2J\x1b[H");
    }
    output.push_str(&format!(
        "Iteration: {iteration:4}  Objective: {objective:10.3}  Volume: {volume:.3}  Change: {change:.3}\n"
    ));

    if interactive {
        for row in 0..densities.nrows() {
            for column in 0..densities.ncols() {
                let density = densities[(row, column)];
                output.push_str(if density > 0.75 {
                    "██"
                } else if density > 0.5 {
                    "▒▒"
                } else if density >= 0.25 {
                    "░░"
                } else if density.is_nan() {
                    "OO"
                } else {
                    "  "
                });
            }
            output.push('\n');
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::format_progress;
    use nalgebra::DMatrix;

    #[test]
    fn redirected_progress_is_one_ascii_line_per_iteration() {
        let densities = DMatrix::from_element(2, 3, 1.0);
        let first = format_progress(1, 123.4567, 0.5, 0.1234, &densities, false);
        let second = format_progress(42, 9.0, 0.5, 0.001, &densities, false);

        assert_eq!(
            first,
            "Iteration:    1  Objective:    123.457  Volume: 0.500  Change: 0.123\n"
        );
        for line in [&first, &second] {
            assert!(line.is_ascii());
            assert!(!line.contains('\x1b'));
            assert!(!line.contains('\t'));
            assert_eq!(line.lines().count(), 1);
            assert!(line.ends_with('\n'));
        }
        for label in ["Objective:", "Volume:", "Change:"] {
            assert_eq!(first.find(label), second.find(label));
        }
    }

    #[test]
    fn interactive_frame_preserves_shading_and_row_order_without_terminal_reset() {
        let densities = DMatrix::from_row_slice(2, 3, &[0.0, 0.25, 0.5, 0.75, 1.0, f64::NAN]);
        let frame = format_progress(1, 10.0, 0.5, 0.2, &densities, true);
        let plain = format_progress(1, 10.0, 0.5, 0.2, &densities, false);

        assert!(frame.starts_with("\x1b[2J\x1b[H"));
        assert!(!frame.contains("\x1bc"));
        assert!(frame.contains(&plain));
        assert!(frame.ends_with("  ░░░░\n▒▒██OO\n"));
        assert_eq!(frame.lines().count(), 3);
    }
}
