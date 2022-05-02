#![warn(clippy::all)]
#![warn(missing_docs)]
#![warn(rustdoc::missing_doc_code_examples)]
#![warn(clippy::missing_docs_in_private_items)]
#![doc = include_str!("../README.md")]

use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::{csc::CscMatrix, factorization::CscCholesky};
mod utils;
use utils::{max, min};

/// The topology optimization function. It takes inputs for the number of elements in the *x*
/// direction (`nelx`), the number of elements in the *y* direction (`nely`), the volume fraction
/// of material to be optimized (`volfrac`), the penalty weighting (`penalty`), and filter
/// radius (`rmin`). It returns a matrix containing the optimized volume of material contained in each
/// cell.
/// ```
///  let x = topopt::top(30, 10, 0.5, 3.0, 1.5, None, None, None, None);
/// ```
pub fn top(
    nelx: usize,
    nely: usize,
    volfrac: f64,
    penalty: f64,
    rmin: f64,
    loads: Option<(DMatrix<f64>, DMatrix<f64>)>,
    boundary: Option<(DMatrix<f64>, DMatrix<f64>)>,
    passive: Option<DMatrix<bool>>,
    active: Option<DMatrix<bool>>,
) -> DMatrix<f64> {
    // INITIALIZE
    let mut x: DMatrix<f64> = DMatrix::from_element(nely, nelx, volfrac);
    let mut xold: DMatrix<f64> = DMatrix::from_element(nely, nelx, volfrac);
    let mut dc: DMatrix<f64> = DMatrix::from_element(nely, nelx, 1.0);
    let mut iter: usize = 0;
    let mut change: f64 = 1.0;
    let mut vol: f64 = 0.0;
    while change > 0.01 {
        iter += 1;
        xold = x.clone();
        // FE-ANALYSIS
        let U = FE(nelx, nely, &x, penalty);

        let mut c = 0.0;
        for ely in 1..=nely {
            for elx in 1..=nelx {
                let n1 = (nely + 1) * (elx - 1) + ely;
                let n2 = (nely + 1) * elx + ely;
                let Ueidx = [
                    2 * n1 - 1,
                    2 * n1,
                    2 * n2 - 1,
                    2 * n2,
                    2 * n2 + 1,
                    2 * n2 + 2,
                    2 * n1 + 1,
                    2 * n1 + 2,
                ];
                let Ue: DVector<f64> = DVector::from_fn(8, |idx, _jdx| U[Ueidx[idx] - 1]);
                let UKEU = (Ue.transpose() * lk() * Ue)[(0, 0)];
                c += x[(ely - 1, elx - 1)].powf(penalty) * UKEU;
                dc[(ely - 1, elx - 1)] =
                    -penalty * x[(ely - 1, elx - 1)].powf(penalty - 1.0) * UKEU;
            }
        }

        // FILTERING OF SENSITIVITIES
        dc = check(nelx, nely, rmin, &x, &dc);

        // % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD

        x = optimality_criteria_update(nelx, nely, &x, volfrac, &dc, &active, &passive);

        // % PRINT RESULTS
        change = (&x - xold).abs().max();
        vol = x.sum() / ((nelx * nely) as f64);

        print!("{esc}c", esc = 27 as char);
        println!("Iter: {iter:04}\tObj: {c:4.3}\tVol: {vol:1.3}\tΔ: {change:1.3}");
        // Print
        for ey in 0..nely {
            for ex in 0..nelx {
                if x[(ey, ex)] > 0.75 {
                    print!("██");
                } else if x[(ey, ex)] > 0.5 {
                    print!("▒▒");
                } else if x[(ey, ex)] >= 0.25 {
                    print!("░░");
                } else if x[(ey, ex)].is_nan() {
                    print!("OO");
                } else {
                    print!("  ");
                }
            }
            print!("\n");
        }
    }
    x
}

#[cfg(test)]
mod top_tests {
    use nalgebra::DMatrix;

    #[test]
    fn test_top() {
        let sol = crate::top(2, 2, 0.5, 3.0, 1.5, None, None, None, None);
        assert!(
            (sol - DMatrix::from_fn(2, 2, |idx, jdx| {
                let x = [[0.5517, 0.3960], [0.5213, 0.5310]];
                x[idx][jdx]
            }))
            .abs()
            .max()
                < 0.001
        )
    }
}

/// Optimality criteria update
pub(crate) fn optimality_criteria_update(
    nelx: usize,
    nely: usize,
    x: &DMatrix<f64>,
    volfrac: f64,
    dc: &DMatrix<f64>,
    active: &Option<DMatrix<bool>>,
    passive: &Option<DMatrix<bool>>,
) -> DMatrix<f64> {
    let mut l1: f64 = 0.0;
    let mut l2: f64 = 100_000.0;
    let delta: f64 = 0.2;
    let mut xnew: DMatrix<f64> = DMatrix::from_element(nelx, nely, 0.0);
    while l2 - l1 > 1e-4 {
        let lmid = 0.5 * (l2 + l1);
        xnew = x.zip_map(dc, |xel, dcel| {
            max(
                0.001,
                max(
                    xel - delta,
                    min(1.0, min(xel + delta, xel * (-dcel / lmid).sqrt())),
                ),
            )
        });

        // Handle active elements
        if let Some(m) = active {
            xnew.zip_apply(m, |mut xel, ael| {
                if ael {
                    *xel = 1.0;
                }
            });
        }

        // Handle passive elements
        if let Some(m) = passive {
            xnew.zip_apply(m, |mut xel, pel| {
                if pel {
                    *xel = 0.001;
                }
            });
        }

        if xnew.sum() - volfrac * (nelx as f64) * (nely as f64) > 0.0 {
            l1 = lmid;
        } else {
            l2 = lmid;
        }
    }
    xnew
}

#[cfg(test)]
mod oc_tests {
    use nalgebra::DMatrix;

    #[test]
    fn test_oc() {
        let oc = crate::optimality_criteria_update(
            3,
            3,
            &DMatrix::from_element(3, 3, 1.0),
            0.5,
            &-DMatrix::from_element(3, 3, 0.1),
            &None,
            &None,
        );
        assert!((oc - DMatrix::from_element(3, 3, 0.8)).max() < 0.001)
    }
}

/// Mesh-independency filter
pub(crate) fn check(
    nelx: usize,
    nely: usize,
    rmin: f64,
    x: &DMatrix<f64>,
    dc: &DMatrix<f64>,
) -> DMatrix<f64> {
    let mut dcn: DMatrix<f64> = DMatrix::from_element(nely, nelx, 0.0);
    for idx in 1..=nelx {
        for jdx in 1..=nely {
            let mut sum = 0.0;
            for kdx in max(idx - rmin.floor() as usize, 1)..=min(idx + rmin.floor() as usize, nelx)
            {
                for ldx in
                    max(jdx - rmin.floor() as usize, 1)..=min(jdx + rmin.floor() as usize, nely)
                {
                    let fac = rmin
                        - (((idx as i32 - kdx as i32).pow(2) + (jdx as i32 - ldx as i32).pow(2))
                            as f64)
                            .sqrt();
                    sum += max(0.0, fac as f64);
                    dcn[(jdx - 1, idx - 1)] +=
                        max(0.0, fac as f64) * x[(ldx - 1, kdx - 1)] * dc[(ldx - 1, kdx - 1)];
                }
            }
            dcn[(jdx - 1, idx - 1)] /= x[(jdx - 1, idx - 1)] * sum;
        }
    }
    dcn
}

#[cfg(test)]
mod check_tests {
    use nalgebra::DMatrix;

    #[test]
    fn test_check() {
        assert!(
            (crate::check(
                3,
                3,
                1.5,
                &DMatrix::from_element(3, 3, 1.0),
                &DMatrix::from_element(3, 3, 0.1)
            ) - DMatrix::from_element(3, 3, 0.1))
            .max()
                < 0.001
        )
    }

    #[test]
    fn test_check_realistic() {
        assert!(
            (crate::check(
                2,
                2,
                1.5,
                &DMatrix::from_element(2, 2, 0.5),
                &DMatrix::from_fn(2, 2, |idx, jdx| {
                    let x = [[-165.5927, -25.7591], [-88.0966, -132.6057]];
                    x[idx][jdx]
                }),
            ) - DMatrix::from_fn(2, 2, |idx, jdx| {
                let x = [[-122.4744, -75.5265], [-109.6200, -104.4333]];
                x[idx][jdx]
            }))
            .abs()
            .max()
                < 0.001
        )
    }
}

/// FE Analysis
pub(crate) fn FE(nelx: usize, nely: usize, x: &DMatrix<f64>, penalty: f64) -> DVector<f64> {
    let KE = lk();
    let mut K: DMatrix<f64> = DMatrix::from_element(
        2 * (nelx + 1) * (nely + 1),
        2 * (nelx + 1) * (nely + 1),
        0.0,
    );
    let mut F: DVector<f64> = DVector::from_element(2 * (nely + 1) * (nelx + 1), 0.0);
    for elx in 1..=nelx {
        for ely in 1..=nely {
            let n1 = (nely + 1) * (elx - 1) + ely;
            let n2 = (nely + 1) * elx + ely;
            let edof = [
                2 * n1 - 1,
                2 * n1,
                2 * n2 - 1,
                2 * n2,
                2 * n2 + 1,
                2 * n2 + 2,
                2 * n1 + 1,
                2 * n1 + 2,
            ];
            for (i, xidx) in edof.iter().enumerate() {
                for (j, yidx) in edof.iter().enumerate() {
                    K[(*xidx - 1, *yidx - 1)] += x[(ely - 1, elx - 1)].powf(penalty) * KE[(i, j)]
                }
            }
        }
    }
    // DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F[1] = -1.0;
    let mut fixeddofs: Vec<usize> = (1..(2 * (nely + 1))).step_by(2).collect();
    fixeddofs.push(2 * (nelx + 1) * (nely + 1));
    fixeddofs.sort_by(|a, b| b.cmp(a));
    for idx in fixeddofs.to_owned() {
        F = F.remove_row(idx - 1);
        K = K.remove_column(idx - 1);
        K = K.remove_row(idx - 1);
    }

    let K_sparse = CscMatrix::from(&K);
    let mut U_as_matrix = CscCholesky::factor(&K_sparse).unwrap().solve(&F);
    let mut U: DVector<f64> =
        DVector::from_fn(U_as_matrix.shape().0, |idx, jdx| U_as_matrix[(idx, 0)]);

    fixeddofs.reverse();
    for idx in fixeddofs.to_owned() {
        U = U.insert_row(idx - 1, 0.0);
    }
    U
}

#[cfg(test)]
mod fe_tests {
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn test_fe() {
        let u_from_matlab: DVector<f64> = DVector::from_vec(vec![
            0.0, -5.6222, 0.0, -4.6222, 1.7222, -1.0000, -2.3222, 0.0,
        ]);

        assert!(
            (crate::FE(1, 1, &DMatrix::from_element(1, 1, 1.0), 10.0) - u_from_matlab)
                .abs()
                .max()
                < 0.001
        );
    }
}

/// Element stiffness matrix
pub(crate) fn lk() -> DMatrix<f64> {
    let elastic_modulus: f64 = 1.0;
    let nu: f64 = 0.3;
    let k = [
        1. / 2. - nu / 6.,
        1. / 8. + nu / 8.,
        -1. / 4. - nu / 12.,
        -1. / 8. + 3. * nu / 8.,
        -1. / 4. + nu / 12.,
        -1. / 8. - nu / 8.,
        nu / 6.,
        1. / 8. - 3. * nu / 8.,
    ];
    DMatrix::from_fn(8, 8, |i, j| {
        let idx = [
            [1, 2, 3, 4, 5, 6, 7, 8],
            [2, 1, 8, 7, 6, 5, 4, 3],
            [3, 8, 1, 6, 7, 4, 5, 2],
            [4, 7, 6, 1, 8, 3, 2, 5],
            [5, 6, 7, 8, 1, 2, 3, 4],
            [6, 5, 4, 3, 2, 1, 8, 7],
            [7, 4, 5, 2, 3, 8, 1, 6],
            [8, 3, 2, 5, 4, 7, 6, 1],
        ];
        elastic_modulus / (1.0 - nu * nu) * k[idx[i][j] - 1]
    })
}

#[cfg(test)]
mod lk_tests {
    use nalgebra::DMatrix;

    #[test]
    fn test_lk() {
        let ke_from_matlab: DMatrix<f64> = DMatrix::from_fn(8, 8, |idx, jdx| {
            let ke = [
                [
                    0.4945, 0.1786, -0.3022, -0.0137, -0.2473, -0.1786, 0.0549, 0.0137,
                ],
                [
                    0.1786, 0.4945, 0.0137, 0.0549, -0.1786, -0.2473, -0.0137, -0.3022,
                ],
                [
                    -0.3022, 0.0137, 0.4945, -0.1786, 0.0549, -0.0137, -0.2473, 0.1786,
                ],
                [
                    -0.0137, 0.0549, -0.1786, 0.4945, 0.0137, -0.3022, 0.1786, -0.2473,
                ],
                [
                    -0.2473, -0.1786, 0.0549, 0.0137, 0.4945, 0.1786, -0.3022, -0.0137,
                ],
                [
                    -0.1786, -0.2473, -0.0137, -0.3022, 0.1786, 0.4945, 0.0137, 0.0549,
                ],
                [
                    0.0549, -0.0137, -0.2473, 0.1786, -0.3022, 0.0137, 0.4945, -0.1786,
                ],
                [
                    0.0137, -0.3022, 0.1786, -0.2473, -0.0137, 0.0549, -0.1786, 0.4945,
                ],
            ];
            ke[idx][jdx]
        });
        assert!((ke_from_matlab - crate::lk()).abs().max() < 0.001);
    }
}
