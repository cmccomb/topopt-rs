/// Minimum function for anything implementing PartialOrd
pub(crate) fn min<T: PartialOrd>(x: T, y: T) -> T {
    if x < y {
        x
    } else {
        y
    }
}

/// Maximum function for anything implementing PartialOrd
pub(crate) fn max<T: PartialOrd>(x: T, y: T) -> T {
    if x > y {
        x
    } else {
        y
    }
}

#[cfg(test)]
mod min_max_tests {

    #[test]
    fn test_min() {
        assert_eq!(crate::min(2, 4), 2);
    }
    #[test]
    fn test_max() {
        assert_eq!(crate::max(2, 4), 4);
    }
}
