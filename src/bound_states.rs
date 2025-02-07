/*

BOUND STATES MODULE

Contains: Calculation of bound states using Numerov method

Last revision: 02/12/2024 by JZ

*/
use nalgebra::{DMatrix, DVector};
use crate::core::{console, HBAR};
use crate::math::{Grid, Integrable};

pub fn calculate_derivative <F: Fn(f64) -> f64> (
    grid:           &Grid,
    potential:      &F,
    e:              f64,
    l:              usize,
    mass:           f64,
) -> (f64,DVector<f64>) {

    let k2: Vec<f64> = grid.points.iter().map(|&r|{
        let k2_i = 2.0 * mass * e / HBAR.powi(2)
            - (l * (l+1)) as f64 / r.powi(2)
            - 2.0 * mass * potential(r) / HBAR.powi(2);
        k2_i
    }).collect();

    let mut psi_r: Vec<f64> = Vec::with_capacity(grid.n);
    psi_r.push(0.0);
    psi_r.push(grid.d);
    psi_r.push(2.0 * (1.0 - 5.0 * grid.d.powi(2) / 12.0 * k2[1]) * grid.d / (1.0 + grid.d.powi(2) / 12.0 * k2[2]));
    for k in 2..grid.n-1 {
        let term2 = 1.0 + grid.d.powi(2) / 12.0 * k2[k + 1];
        let term1 = 1.0 - 5.0 * grid.d.powi(2) / 12.0 * k2[k];
        let term0 = 1.0 + grid.d.powi(2) / 12.0 * k2[k - 1];
        psi_r.push((2.0 * term1 * psi_r[k] - term0 * psi_r[k-1]) / term2);

        if psi_r[k+1] / psi_r[k] < 1.0 {
            break;
        }
    }
    psi_r.shrink_to_fit();
    let i_break = psi_r.len() - 1;
    let norm_r = 1.0/psi_r[i_break-1];
    psi_r.iter_mut().for_each(|x| *x *= norm_r);;
    let mut psi_i= vec![0.0;grid.n-psi_r.len()+3];
    psi_i[0] = (-(-k2[grid.n-1]).sqrt()*grid.points[grid.n-1]).exp();
    psi_i[1] = (-(-k2[grid.n-2]).sqrt()*grid.points[grid.n-2]).exp();
    let mut index = 0;
    loop {
        if index + 4 >= psi_i.len() || psi_i[1 + index].abs() > 0.0 {
            break;
        }
        psi_i[2 + index] = (-(-k2[grid.n-3-index]).sqrt()*grid.points[grid.n-3-index]).exp();
        psi_i[3 + index] = (-(-k2[grid.n-4-index]).sqrt()*grid.points[grid.n-4-index]).exp();

        index += 2;
    }
    for k in 1+index..(grid.n-i_break+1) {
        let term2 = 1.0 + grid.d.powi(2) / 12.0 * k2[grid.n-2-k];
        let term1 = 1.0 - 5.0 * grid.d.powi(2) / 12.0 * k2[grid.n-1-k];
        let term0 = 1.0 + grid.d.powi(2) / 12.0 * k2[grid.n-k];
        psi_i[k+1] = ((2.0 * term1 * psi_i[k] - term0 * psi_i[k-1]) / term2);
    }
    let norm_i = 1.0/psi_i[grid.n-i_break];
    psi_i.iter_mut().for_each(|x| *x *= norm_i);
    psi_i.reverse();

    if psi_r[i_break-1] != psi_i[1] {
        console(&format!("[ERROR]: psi_r and psi_i are not correctly normalized at energy {} au.",&e));
    }

    let derivative_diff = (psi_i[2]-psi_i[0] - psi_r[i_break]+psi_r[i_break-2])/(2.0*grid.d);
    psi_i.drain(0..3);
    psi_r.extend(psi_i);
    let psi = DVector::from_vec(psi_r);
    (derivative_diff, psi)
}


pub fn calculate_bound_states(
    grid:           &Grid,
    potential:      impl Fn(f64) -> f64,
    l:              usize,
    mass:           f64,
    ne:             usize,
) -> (DVector<f64>,DMatrix<f64>,DVector<f64>,DVector<f64>) {

    let effective_potential: Vec<f64> = grid.points
        .iter()
        .map(|&r| potential(r) + (l * (l + 1)) as f64 / (2.0 * mass * r.powi(2)))
        .collect();

    let emin = effective_potential
        .iter()
        .filter(|&&x| x.is_finite())
        .copied()
        .reduce(f64::min)
        .unwrap_or(-f64::INFINITY);

    let emax = effective_potential[effective_potential.len() - 1];
    console(&format!("Searching for bound states in energy range [{emin},{emax}] au"));

    let energies = Grid::new_n(emin,emax,ne);
    let derivative_values: Vec<f64> = energies.points.iter().map(|&e| {
        let (derivative, _psi) = calculate_derivative(
            &grid,
            &potential,
            e,
            l,
            mass
        );
        derivative
    }).collect();

    let mut bound_energies: Vec<f64> = Vec::new();

    for (i, window) in derivative_values.windows(2).enumerate() {
        if window[0] * window[1] <= 0.0 {
            let mut e_left = energies.points[i];
            let mut e_right = energies.points[i + 1];
            let mut converged = false;

            for _ in 0..1000 {
                let e_mid = (e_left + e_right) / 2.0;
                let (derivative_mid, _) = calculate_derivative(
                    &grid,
                    &potential,
                    e_mid,
                    l,
                    mass);

                if derivative_mid.abs() <= 2.0 * grid.d.powi(2) {
                    bound_energies.push(e_mid);
                    converged = true;
                    break;
                }

                if derivative_mid * window[0] <= 0.0 {
                    e_right = e_mid;
                } else {
                    e_left = e_mid;
                }
            }

            if !converged {
                console(&format!("[ERROR]: Bisection did not converge at energy {} au", e_left));
            }
        }
    }
    let final_energies = DVector::from_vec(bound_energies);
    let mut eigenfunctions = DMatrix::zeros(final_energies.len(),grid.n);
    for (i,&e) in final_energies.iter().enumerate() {
        let (_dif, psi) = calculate_derivative(
            &grid,
            &potential,
            e,
            l,
            mass
        );
        let norm = psi.square_integrate(grid.d);
        let normalized_psi =  psi * (1.0 / norm.sqrt());
        eigenfunctions.row_mut(i).copy_from_slice(&normalized_psi.as_slice());
    }
    let matching_function = DVector::from_vec(derivative_values);
    let computed_energies = DVector::from_vec(energies.points);
    (final_energies,eigenfunctions,computed_energies,matching_function)
}