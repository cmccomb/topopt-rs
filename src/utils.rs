use nalgebra::DMatrix;

/// Elementwise multiplication
pub(crate) fn mult(x: DMatrix<f64>, y: DMatrix<f64>) -> DMatrix<f64> {
    x.zip_map(&y, |x, y| x * y)
}

/// Elementwise minimum
pub(crate) fn min(x: DMatrix<f64>, y: DMatrix<f64>) -> DMatrix<f64> {
    x.zip_map(&y, |x, y| tmin(x, y))
}

/// Elementwise maximum
pub(crate) fn max(x: DMatrix<f64>, y: DMatrix<f64>) -> DMatrix<f64> {
    -min(-x, -y)
}

/// Minimum function for anything implementing PartialOrd
pub(crate) fn tmin<T: PartialOrd>(x: T, y: T) -> T {
    if x < y {
        x
    } else {
        y
    }
}

/// Maximum function for anything implementing PartialOrd
pub(crate) fn tmax<T: PartialOrd>(x: T, y: T) -> T {
    if x > y {
        x
    } else {
        y
    }
}

pub(crate) fn nancount(matrix: &DMatrix<f64>) -> usize {
    let mut counter = 0;
    matrix.map(|x| {
        if x.is_nan() {
            counter += 1
        }
    });
    counter
}

#[cfg(test)]
mod min_max_tests {
    use nalgebra::DMatrix;

    #[test]
    fn test_tmin() {
        assert_eq!(crate::tmin(2, 4), 2);
    }
    #[test]
    fn test_tmax() {
        assert_eq!(crate::tmax(2, 4), 4);
    }
    #[test]
    fn test_min() {
        let m1 = DMatrix::from_vec(2, 2, vec![1., 2., 1., 2.]);
        let m2 = DMatrix::from_vec(2, 2, vec![2., 1., 2., 1.]);
        assert_eq!(crate::min(m1, m2), DMatrix::from_element(2, 2, 1.0));
    }
    #[test]
    fn test_max() {
        let m1 = DMatrix::from_vec(2, 2, vec![1., 2., 1., 2.]);
        let m2 = DMatrix::from_vec(2, 2, vec![2., 1., 2., 1.]);
        assert_eq!(crate::max(m1, m2), DMatrix::from_element(2, 2, 2.0));
    }
}
