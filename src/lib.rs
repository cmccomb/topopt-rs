#![warn(clippy::all)]
#![warn(missing_docs)]
#![warn(rustdoc::missing_doc_code_examples)]
#![warn(clippy::missing_docs_in_private_items)]
#![doc = include_str!("../README.md")]

use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::{csc::CscMatrix, factorization::CscCholesky};
mod utils;
use utils::{max, min};

/// The topology optimization solver.
///
/// It takes the following inputs
/// - `nelx`, the number of elements in the *x* direction
/// - `nely`, the number of elements in the *x* direction
/// - `volfrac`, the volume fraction of material to be optimized for
/// - `penalty` the penalty weight
/// - `rmin`, the filter radius
/// - `loads`, an optional `nelx+1`-by-`nely+1` matrix of `f64` tuples indicating the load applied to each node in *x*/*y* pairs.
/// - `boundary`, an optional `nelx+1`-by-`nely+1` matrix of `bool` tuples indicating the degrees of freedom of each node in *x*/*y* pairs
/// - `passive`, an optional `nelx`-by-`nely` matrix of `bool`s indicating elements which should be active (always void)
/// - `active`, an optional `nelx`-by-`nely` matrix of `bool`s indicating elements which should be active (always filled)
///
/// It returns a matrix of size `nelx`-by-`nely` containing the optimized volume of material contained in each
/// cell.
/// ```
///  let x = topopt::top(30, 10, 0.5, 3.0, 1.5, None, None, None, None);
/// ```
pub fn top(
    nelx: usize, // asdfasdf
    nely: usize,
    volfrac: f64,
    penalty: f64,
    rmin: f64,
    loads: Option<DMatrix<(f64, f64)>>,
    boundary: Option<DMatrix<(bool, bool)>>,
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
        let U = finite_element(nelx, nely, &x, penalty, &loads, &boundary);

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

pub fn solve(settings: Settings) -> DMatrix<f64> {
    top(
        settings.nelx,
        settings.nely,
        settings.volume_fraction,
        settings.penalty_weight,
        settings.filter_radius,
        Some(settings.loads),
        Some(settings.boundary),
        Some(settings.passive),
        Some(settings.active),
    )
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
pub(crate) fn finite_element(
    nelx: usize,
    nely: usize,
    x: &DMatrix<f64>,
    penalty: f64,
    loads: &Option<DMatrix<(f64, f64)>>,
    boundary: &Option<DMatrix<(bool, bool)>>,
) -> DVector<f64> {
    let KE = lk();
    let mut K: DMatrix<f64> = DMatrix::from_element(
        2 * (nelx + 1) * (nely + 1),
        2 * (nelx + 1) * (nely + 1),
        0.0,
    );
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
    let mut F: DVector<f64> = DVector::from_element(2 * (nely + 1) * (nelx + 1), 0.0);
    match loads {
        None => {
            F[1] = -1.0;
        }
        Some(load_matrix) => {
            let mut counter: usize = 0;
            for idx in 0..nelx {
                for jdx in 0..nely {
                    F[counter] = load_matrix[(idx, jdx)].0;
                    F[counter + 1] = load_matrix[(idx, jdx)].1;
                    counter += 2;
                }
            }
        }
    }

    let mut fixeddofs: Vec<usize> = vec![];
    match boundary {
        None => {
            fixeddofs = (1..(2 * (nely + 1))).step_by(2).collect();
            fixeddofs.push(2 * (nelx + 1) * (nely + 1));
        }
        Some(boundary_matrix) => {
            let mut counter: usize = 1;
            for idx in 0..=nelx {
                for jdx in 0..=nely {
                    if boundary_matrix[(idx, jdx)].0 {
                        fixeddofs.push(counter);
                    }

                    if boundary_matrix[(idx, jdx)].1 {
                        fixeddofs.push(counter + 1);
                    }

                    counter += 2;
                }
            }
        }
    }
    fixeddofs.sort_by(|a, b| b.cmp(a));

    // Do magic
    for idx in fixeddofs.to_owned() {
        F = F.remove_row(idx - 1);
        K = K.remove_column(idx - 1);
        K = K.remove_row(idx - 1);
    }

    // Solve matrix
    let K_sparse = CscMatrix::from(&K);
    let mut U_as_matrix = CscCholesky::factor(&K_sparse).unwrap().solve(&F);
    let mut U: DVector<f64> =
        DVector::from_fn(U_as_matrix.shape().0, |idx, jdx| U_as_matrix[(idx, 0)]);

    // Undo magic
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
    fn test_finite_element() {
        let u_from_matlab: DVector<f64> = DVector::from_vec(vec![
            0.0, -5.6222, 0.0, -4.6222, 1.7222, -1.0000, -2.3222, 0.0,
        ]);

        assert!(
            (crate::finite_element(1, 1, &DMatrix::from_element(1, 1, 1.0), 10.0, &None, &None)
                - u_from_matlab)
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

///
#[derive(Clone)]
pub struct Settings {
    nelx: usize,
    nely: usize,
    volume_fraction: f64,
    filter_radius: f64,
    penalty_weight: f64,
    loads: DMatrix<(f64, f64)>,
    boundary: DMatrix<(bool, bool)>,
    passive: DMatrix<bool>,
    active: DMatrix<bool>,
}

impl Default for Settings {
    /// ```
    /// use topopt::Settings;
    /// Settings::default();
    /// ```
    fn default() -> Self {
        Self {
            nelx: 60,
            nely: 20,
            volume_fraction: 0.5,
            filter_radius: 1.5,
            penalty_weight: 3.0,
            loads: DMatrix::from_fn(61, 21, |idx, jdx| {
                if idx == 0 && jdx == 0 {
                    (0.0, -1.0)
                } else {
                    (0.0, 0.0)
                }
            }),
            boundary: DMatrix::from_fn(61, 21, |idx, jdx| {
                if idx == 0 {
                    (true, false)
                } else if idx == 60 && jdx == 20 {
                    (false, true)
                } else {
                    (false, false)
                }
            }),
            passive: DMatrix::from_element(20, 60, false),
            active: DMatrix::from_element(20, 60, false),
        }
    }
}

impl Settings {
    /// ```
    /// use topopt::Settings;
    /// Settings::new(60, 20, 0.5);
    /// ```
    pub fn new(nelx: usize, nely: usize, volume_fraction: f64) -> Self {
        Self {
            nelx,
            nely,
            volume_fraction,
            filter_radius: 1.5,
            penalty_weight: 3.0,
            loads: DMatrix::from_element(nelx + 1, nely + 1, (0.0, 0.0)),
            boundary: DMatrix::from_element(nelx + 1, nely + 1, (false, false)),
            passive: DMatrix::from_element(nely, nelx, false),
            active: DMatrix::from_element(nely, nelx, false),
        }
    }

    pub fn with_filter_radius(&mut self, filter_radius: f64) -> Self {
        self.filter_radius = filter_radius;
        self.clone()
    }

    pub fn with_penalty_weight(&mut self, penalty_weight: f64) -> Self {
        self.penalty_weight = penalty_weight;
        self.clone()
    }

    pub fn with_active_elements(&mut self, mask: DMatrix<bool>) -> Self {
        self.active = mask;
        self.clone()
    }
    pub fn with_passive_elements(&mut self, mask: DMatrix<bool>) -> Self {
        self.passive = mask;
        self.clone()
    }
}

/// # Modifying the Boundary Conditions
/// This group of methods provides utilities for changing the boundary conditions. These methods are applied like so:
///```
/// # use topopt::Settings;
/// let settings = Settings::new(120, 30, 0.5)
///     .with_bottom_right_bc(false, true)
///     .with_bottom_left_bc(false, true);
/// ```
impl Settings {
    pub fn with_left_bc(&mut self, x: bool, y: bool) -> Self {
        for idx in 0..=self.nelx {
            for jdx in 0..=self.nely {
                if idx == 0 {
                    self.boundary[(idx, jdx)] = (x, y);
                }
            }
        }
        self.clone()
    }

    pub fn with_right_bc(&mut self, x: bool, y: bool) -> Self {
        for idx in 0..=self.nelx {
            for jdx in 0..=self.nely {
                if idx == self.nelx {
                    self.boundary[(idx, jdx)] = (x, y);
                }
            }
        }
        self.clone()
    }

    pub fn with_top_bc(&mut self, x: bool, y: bool) -> Self {
        for idx in 0..=self.nelx {
            for jdx in 0..=self.nely {
                if jdx == 0 {
                    self.boundary[(idx, jdx)] = (x, y);
                }
            }
        }
        self.clone()
    }

    pub fn with_bottom_bc(&mut self, x: bool, y: bool) -> Self {
        for idx in 0..=self.nelx {
            for jdx in 0..=self.nely {
                if jdx == self.nely {
                    self.boundary[(idx, jdx)] = (x, y);
                }
            }
        }
        self.clone()
    }

    pub fn with_vertical_midline_bc(&mut self, x: bool, y: bool) -> Self {
        for idx in 0..=self.nelx {
            for jdx in 0..=self.nely {
                if idx == self.nelx / 2 {
                    self.boundary[(idx, jdx)] = (x, y);
                }
            }
        }
        self.clone()
    }

    pub fn with_horizontal_midline_bc(&mut self, x: bool, y: bool) -> Self {
        for idx in 0..=self.nelx {
            for jdx in 0..=self.nely {
                if jdx == self.nely / 2 {
                    self.boundary[(idx, jdx)] = (x, y);
                }
            }
        }
        self.clone()
    }

    pub fn with_bottom_right_bc(&mut self, x: bool, y: bool) -> Self {
        self.boundary[(self.nelx, self.nely)] = (x, y);
        self.clone()
    }

    pub fn with_bottom_left_bc(&mut self, x: bool, y: bool) -> Self {
        self.boundary[(0, self.nely)] = (x, y);
        self.clone()
    }

    pub fn with_top_right_bc(&mut self, x: bool, y: bool) -> Self {
        self.boundary[(self.nelx, 0)] = (x, y);
        self.clone()
    }

    pub fn with_top_left_bc(&mut self, x: bool, y: bool) -> Self {
        self.boundary[(0, 0)] = (x, y);
        self.clone()
    }

    pub fn with_bc(&mut self, boundary: DMatrix<(bool, bool)>) -> Self {
        self.boundary = boundary;
        self.clone()
    }

    pub fn set_bc(&mut self, idx: usize, jdx: usize, x: bool, y: bool) -> Self {
        self.boundary[(idx, jdx)] = (x, y);
        self.clone()
    }
}
/// # Modifying the Load Case
/// This group of methods provides utilities for changing the load case. These methods are applied like so:
///```
/// # use topopt::Settings;
/// let settings = Settings::new(120, 30, 0.5)
///     .with_bottom_right_bc(false, true)
///     .with_bottom_left_bc(false, true)
///     .with_top_middle_load(0.0, -1.0);
/// ```
impl Settings {
    pub fn with_bottom_right_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(self.nelx, self.nely)] = (x, y);
        self.clone()
    }

    pub fn with_bottom_left_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(0, self.nely)] = (x, y);
        self.clone()
    }

    pub fn with_top_right_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(self.nelx, 0)] = (x, y);
        self.clone()
    }

    pub fn with_top_left_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(0, 0)] = (x, y);
        self.clone()
    }

    pub fn with_top_middle_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(self.nelx / 2, 0)] = (x, y);
        self.clone()
    }

    pub fn with_bottom_middle_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(self.nelx / 2, self.nely)] = (x, y);
        self.clone()
    }

    pub fn with_right_middle_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(self.nelx, self.nely / 2)] = (x, y);
        self.clone()
    }

    pub fn with_left_middle_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(0, self.nely / 2)] = (x, y);
        self.clone()
    }

    pub fn with_centered_load(&mut self, x: f64, y: f64) -> Self {
        self.loads[(self.nelx / 2, self.nely / 2)] = (x, y);
        self.clone()
    }

    pub fn with_loads(&mut self, loads: DMatrix<(f64, f64)>) -> Self {
        self.loads = loads;
        self.clone()
    }

    pub fn set_load(&mut self, idx: usize, jdx: usize, x: bool, y: bool) -> Self {
        self.boundary[(idx, jdx)] = (x, y);
        self.clone()
    }
}
