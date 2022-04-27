use nalgebra::{DMatrix, DVector};
mod utils;
use utils::{max, min, mult, tmax, tmin};

/// ```
///  let x = topopt::top(30, 10, 0.5, 3.0, 1.5);
/// ```
pub fn top(nelx: usize, nely: usize, volfrac: f32, penalty: f32, rmin: f32) -> DMatrix<f32> {
    // INITIALIZE
    let mut x: DMatrix<f32> = DMatrix::from_element(nely, nelx, volfrac);
    let mut xold: DMatrix<f32> = DMatrix::from_element(nely, nelx, volfrac);
    let mut dc: DMatrix<f32> = DMatrix::from_element(nely, nelx, 1.0);
    let mut iter: usize = 0;
    let mut change: f32 = 1.0;
    // START ITERATION
    while change > 0.01 {
        iter += 1;
        xold = x.clone();
        // FE-ANALYSIS
        let U = FE(nelx, nely, x.clone(), penalty);
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
                let Ue: DVector<f32> = DVector::from_fn(8, |idx, _jdx| U[Ueidx[idx] - 1]);
                let UKEU = (Ue.transpose() * lk() * Ue)[(0, 0)];
                c += x[(ely - 1, elx - 1)].powf(penalty) * UKEU;
                dc[(ely - 1, elx - 1)] =
                    -penalty * x[(ely - 1, elx - 1)].powf(penalty - 1.0) * UKEU;
            }
        }

        // FILTERING OF SENSITIVITIES
        dc = check(nelx, nely, rmin, x.clone(), dc.clone());

        // % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
        x = optimality_criteria_update(nelx, nely, x.clone(), volfrac, dc.clone());
        // % PRINT RESULTS
        change = (x.clone() - xold).abs().max();
        let vol = x.sum() / ((nelx * nely) as f32);
        println!("It.: {iter} Obj.: {c} Vol.: {vol} ch.: {change}");
    }
    x
}

#[cfg(test)]
mod top_tests {
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn test_top() {
        let sol = crate::top(2, 2, 0.5, 3.0, 1.5);
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
    x: DMatrix<f32>,
    volfrac: f32,
    dc: DMatrix<f32>,
) -> DMatrix<f32> {
    let mut l1: f32 = 0.0;
    let mut l2: f32 = 100_000.0;
    let delta: f32 = 0.2;
    let mut xnew = x.clone();
    while l2 - l1 > 1e-4 {
        let lmid = 0.5 * (l2 + l1);
        let dc_transform: DMatrix<f32> = dc.scale(-1.0 / lmid).map(|x| x.sqrt());
        xnew = max(
            DMatrix::from_element(nely, nelx, 0.001),
            max(
                x.add_scalar(-delta),
                min(
                    DMatrix::from_element(nely, nelx, 1.0),
                    min(x.add_scalar(delta), mult(x.clone(), dc_transform)),
                ),
            ),
        );
        if xnew.sum() - volfrac * (nelx as f32) * (nely as f32) > 0.0 {
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
            DMatrix::from_element(3, 3, 1.0),
            0.5,
            -DMatrix::from_element(3, 3, 0.1),
        );
        assert!((oc - DMatrix::from_element(3, 3, 0.8)).max() < 0.001)
    }
}

/// Mesh-independency filter
pub(crate) fn check(
    nelx: usize,
    nely: usize,
    rmin: f32,
    x: DMatrix<f32>,
    dc: DMatrix<f32>,
) -> DMatrix<f32> {
    let mut dcn: DMatrix<f32> = DMatrix::from_element(nely, nelx, 0.0);
    for idx in 1..=nelx {
        for jdx in 1..=nely {
            let mut sum = 0.0;
            for kdx in
                tmax(idx - rmin.floor() as usize, 1)..=tmin(idx + rmin.floor() as usize, nelx)
            {
                for ldx in
                    tmax(jdx - rmin.floor() as usize, 1)..=tmin(jdx + rmin.floor() as usize, nely)
                {
                    let fac = rmin
                        - (((idx as i32 - kdx as i32).pow(2) + (jdx as i32 - ldx as i32).pow(2))
                            as f32)
                            .sqrt();
                    sum = sum + tmax(0.0, fac as f32);
                    dcn[(jdx - 1, idx - 1)] +=
                        tmax(0.0, fac as f32) * x[(ldx - 1, kdx - 1)] * dc[(ldx - 1, kdx - 1)];
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
                DMatrix::from_element(3, 3, 1.0),
                DMatrix::from_element(3, 3, 0.1)
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
                DMatrix::from_element(2, 2, 0.5),
                DMatrix::from_fn(2, 2, |idx, jdx| {
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
pub(crate) fn FE(nelx: usize, nely: usize, x: DMatrix<f32>, penalty: f32) -> DVector<f32> {
    let KE = lk();
    let mut K: DMatrix<f32> = DMatrix::from_element(
        2 * (nelx + 1) * (nely + 1),
        2 * (nelx + 1) * (nely + 1),
        0.0,
    );
    let mut F: DVector<f32> = DVector::from_element(2 * (nely + 1) * (nelx + 1), 0.0);
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
    // This is not in keeping with tradition at all
    F[1] = -1.0;
    let mut fixeddofs: Vec<usize> = (1..(2 * (nely + 1))).step_by(2).collect();
    fixeddofs.push(2 * (nelx + 1) * (nely + 1));
    fixeddofs.sort_by(|a, b| b.cmp(a));
    for idx in fixeddofs.to_owned() {
        F = F.remove_row(idx - 1);
        K = K.remove_column(idx - 1);
        K = K.remove_row(idx - 1);
    }

    let mut U: DVector<f32> = K.try_inverse().expect("Cannot invert K") * F;

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
        let u_from_matlab: DVector<f32> = DVector::from_vec(vec![
            0.0, -5.6222, 0.0, -4.6222, 1.7222, -1.0000, -2.3222, 0.0,
        ]);

        assert!(
            (crate::FE(1, 1, DMatrix::from_element(1, 1, 1.0), 10.0) - u_from_matlab)
                .abs()
                .max()
                < 0.001
        );
    }
}

/// Element stiffness matrix
pub(crate) fn lk() -> DMatrix<f32> {
    let elastic_modulus: f32 = 1.0;
    let nu: f32 = 0.3;
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
        let ke_from_matlab: DMatrix<f32> = DMatrix::from_fn(8, 8, |idx, jdx| {
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
