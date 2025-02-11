/*

SCATTERING MODULE

Contains: Scattering in general potential, eigenstates calculation

Last revision: 07/02/2025 by JZ

*/
use nalgebra::{DMatrix, DVector};
use crate::core::{console, HBAR, PI};
use crate::math::{riccati_bessel, riccati_neumann, Grid};

pub fn calculate_phase_shift(
    grid:           &Grid,
    potential:      impl Fn(f64) -> f64,
    energies:       &Grid,
    ls:             &[usize],
    mass:           f64,
) -> DMatrix<f64> {
    let mut phase_shifts = DMatrix::zeros(energies.n,ls.len());
    let shift:usize = 30;

    for (i,&e) in energies.points.iter().enumerate() {
        for (j, &l) in ls.iter().enumerate() {

            let k2: Vec<f64> = grid.points.iter().map(|&r|{
                let k2_i = 2.0 * mass * e / HBAR.powi(2)
                                - (l * (l+1)) as f64 / r.powi(2)
                                - 2.0 * mass * potential(r) / HBAR.powi(2);
                k2_i
            }).collect();

            let mut psi = vec![0.0;grid.n];
            psi[1] = grid.d;
            psi[2] = 2.0 * (1.0 - 5.0 * grid.d.powi(2) / 12.0 * k2[1]) * psi[1] / (1.0 + grid.d.powi(2) / 12.0 * k2[2]);

            for k in 2..(grid.n-1) {
                let term2 = 1.0 + grid.d.powi(2) / 12.0 * k2[k+1];
                let term1 = 1.0 - 5.0 * grid.d.powi(2) / 12.0 * k2[k];
                let term0 = 1.0 + grid.d.powi(2) / 12.0 * k2[k-1];
                psi[k+1] = (2.0 * term1 * psi[k] - term0 * psi[k-1]) / term2;
            }

            let i2 = grid.n - 1;
            let i1 = grid.n - 1 - shift;
            let r2 = grid.points[i2];
            let r1 = grid.points[i1];

            let g = psi[i1]/psi[i2];
            if g.is_nan() {console(&format!("[ERROR]: g is NaN for e = {} a.u. and l = {}",e,l))}
            let k = (2.0 * mass * e).sqrt() / HBAR;

            let delta = ((riccati_bessel(l,k*r1) - g * riccati_bessel(l,k*r2)) / (g * riccati_neumann(l,k*r2) - riccati_neumann(l,k*r1))).atan();
            phase_shifts[(i,j)] = delta;
        }
    }

    phase_shifts
}



pub fn calculate_cross_sections(
    phase_shifts:   &DMatrix<f64>,
    energies:       &Grid,
    ls:             &[usize],
) -> (DMatrix<f64>, DVector<f64>) {
    let mut cross_sections_l = DMatrix::zeros(energies.n,ls.len());
    let mut cross_sections = DVector::zeros(energies.n);

    for (i,&e) in energies.points.iter().enumerate() {
        let normalisation = 4.0 * PI  / (2.0 * e / HBAR.powi(2));

        for (j,&l) in ls.iter().enumerate() {
            let phase_shift = phase_shifts[(i,j)];
            let cross_section_l = normalisation * (2 * l + 1) as f64 * phase_shift.sin().powi(2);
            cross_sections_l[(i,j)] = cross_section_l;
            cross_sections[i] += cross_section_l;
        }
    }
    (cross_sections_l, cross_sections)
}