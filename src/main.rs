use std::fs::File;
use nalgebra::Complex;
use crate::core::{console, Settings, PHYS_H0, PHYS_ME, PHYS_U};
use crate::math::{riccati_bessel, riccati_neumann, Grid};
use std::io::{Write};
use num_traits::{Float, Pow};
use std::result::Result;
use crate::bound_states::{calculate_bound_states, calculate_derivative};
use crate::scattering::{calculate_cross_sections, calculate_phase_shift};

mod core;
mod math;
mod scattering;
mod bound_states;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let settings = Settings::initialize("settings.toml")?;
    let x = Grid::new_n(settings.rmin,settings.rmax,settings.nr);
    let e = Grid::new_d(settings.emin,settings.emax,settings.de);
    let l: Vec<usize> = (0..=settings.lmax).collect();
    let m = 1.0;
    //let m = 14.0 * PHYS_U / PHYS_ME / 2.0;

   /* let v = |r: f64| -> f64 {
        let v0 = 0.75102;
        let alpha = 1.15350;
        let r0 = 2.01943;
        v0 * ((-2.0*alpha*(r-r0)).exp() - 2.0*(-alpha*(r-r0)).exp())
    };

    */

    let v = |r: f64| -> f64 {
        let a = 1.0;
        let alpha = 4.8;
        return if { r <= a } {
            -alpha.powi(2) / (2.0 * m * a.powi(2))
        } else { 0.0 }
    };



    /*
    let v = |r: f64| -> f64 {
        let r0 = 2.0787;
        let a = 2.42;
        let b = -1.07 / PHYS_H0;
        let d = 12.2 / PHYS_H0;

        d * (-2.0 * a * (r - r0) /r0).exp() - 2.0 * d * (-a * (r - r0) /r0).exp() + b
    };
    */

    // Potential output \\
    let mut file = File::create("V.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "R [au]", "V(R) [eV]")?;
    for i in 0..x.n {
        write!(file, "{:20.12E}", x.points[i])?;
        write!(file, "{:20.12E}", v(x.points[i]))?;
        writeln!(file)?;
    }
    console("Data successfully written to V.txt");


    // Phase-shifts output \\
    let mut file = File::create("phase_shifts.txt")?;
    let phase_shifts = calculate_phase_shift(&x,v,&e,&l,m);
    writeln!(file, "{:<3}{:>17}{:>20}{:>20}{:>20}{:>20}", "#", "E [au]", "delta_0 [1]", "delta_1 [1]", "delta_2 [1]", "delta_3 [1]")?;
    for i in 0..e.n {
        write!(file, "{:20.12E}", e.points[i])?;
        for j in 0..l.len() {
            write!(file, "{:20.12E}", phase_shifts[(i,j)])?;
        }
        writeln!(file)?;
    }
    console("Data successfully written to phase_shifts.txt");


    // Cross sections output \\
    let mut file = File::create("cross_sections.txt")?;
    let (cross_sections,cross_section_tot) = calculate_cross_sections(&phase_shifts,&e,&l);
    writeln!(file, "{:<3}{:>17}{:>20}{:>20}{:>20}{:>20}{:>20}", "#", "E [au]", "sigma_0 [au^2]", "sigma_1 [au^2]", "sigma_2 [au^2]", "sigma_3 [au^2]", "sigma_tot [au^2]")?;
    for i in 0..e.n {
        write!(file, "{:20.12E}", e.points[i])?;
        for j in 0..l.len() {
            write!(file, "{:20.12E}", cross_sections[(i,j)])?;
        }
        write!(file, "{:20.12E}", cross_section_tot[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to cross_sections.txt");

    // Bound states calculation \\
    let (energies,functions,calculated_e,matching_f) = calculate_bound_states(&x,v,0,m,settings.bound_ne);

    let mut file = File::create("energies.txt")?;
    writeln!(file, "{:<3}{:>1}{:>20}", "#", "n", "E_n [a.u.]")?;
    for i in 0..energies.len() {
        write!(file, "{:4.0}", i)?;
        write!(file, "{:20.12E}", energies[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to energies.txt");

    let mut file = File::create("matching_function.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "e [au]", "f [a.u.]")?;
    for i in 0..calculated_e.len() {
        write!(file, "{:20.12E}", calculated_e[i])?;
        write!(file, "{:20.12E}", matching_f[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to energies.txt");

    let mut file = File::create("eigenfunctions.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}{:>20}{:>20}{:>20}{:>20}", "#", "R [au]", "psi_0 [a.u.]", "psi_1 [a.u.]", "psi_2 [a.u.]", "psi_3 [a.u.]", "psi_... [a.u.]")?;
    for i in 0..x.n {
        write!(file, "{:20.12E}", x.points[i])?;
        for j in 0..energies.len() {
            write!(file, "{:20.12E}", functions[(j,i)])?;
        }
        writeln!(file)?;
    }
    console("Data successfully written to eigenfunctions.txt");

    Ok(())
}
