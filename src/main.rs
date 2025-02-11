use std::fs;
use std::fs::File;
use nalgebra::{DVector};
use crate::core::{console, Settings};
use crate::math::{Grid};
use std::io::{Write};
use std::result::Result;
use crate::bound_states::{calculate_bound_states};
use crate::scattering::{calculate_cross_sections, calculate_phase_shift};

mod core;
mod math;
mod scattering;
mod bound_states;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    fs::create_dir_all("Test")?;
    fs::create_dir_all("Task 1")?;
    fs::create_dir_all("Task 2")?;






    // Testing output \\

    let settings = Settings::initialize("settings.toml")?;
    let r = Grid::new_n(settings.rmin1,settings.rmax1,settings.nr1);
    let r_bound = Grid::new_n(settings.rmin1,settings.rmax1 * settings.bound_mult as f64,settings.nr1 * settings.bound_mult);
    let e = Grid::new_d(settings.emin1,settings.emax1,settings.de1);
    let e_zero = Grid::new_log(std::f64::EPSILON,settings.de1,settings.zero_limit);
    let ls: Vec<usize> = settings.ls1;
    let m = settings.mass1;

   let test_pot = |_: f64| -> f64 {
       0.0
    };

    // Potential output \\
    let mut file = File::create("Test/test_potential.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "R [au]", "V(R) [au]")?;
    for i in 0..r.n {
        write!(file, "{:20.12E}", r.points[i])?;
        write!(file, "{:20.12E}", test_pot(r.points[i]))?;
        writeln!(file)?;
    }
    console("Data successfully written to Test/test_potential.txt");

    // Phase-shifts output \\
    let mut file = File::create("Test/phase_shifts.txt")?;
    let phase_shifts = calculate_phase_shift(&r,test_pot,&e,&ls,m);
    write!(file, "{:<3}{:>17}", "#", "E [au]")?;
    for l in ls.iter() {
        write!(file,"{:>20}",&format!("delta_{} [1]",l))?;
    }
    writeln!(file);
    for i in 0..e.n {
        write!(file, "{:20.12E}", e.points[i])?;
        for j in 0..ls.len() {
            write!(file, "{:20.12E}", phase_shifts[(i,j)])?;
        }
        writeln!(file)?;
    }
    console("Data successfully written to Test/phase_shifts.txt");

    // Cross sections output \\
    let mut file = File::create("Test/cross_sections.txt")?;
    let (cross_sections,cross_section_tot) = calculate_cross_sections(&phase_shifts,&e,&ls);
    write!(file, "{:<3}{:>17}", "#", "E [au]")?;
    for l in ls.iter() {
        write!(file,"{:>20}",&format!("sigma_{} [au^2]",l))?;
    }
    writeln!(file,"{:>20}","sigma_tot [au^2]")?;
    for i in 0..e.n {
        write!(file, "{:20.12E}", e.points[i])?;
        for j in 0..ls.len() {
            write!(file, "{:20.12E}", cross_sections[(i,j)])?;
        }
        write!(file, "{:20.12E}", cross_section_tot[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to Test/cross_sections.txt");




    // Task 1 output \\


    for alpha in settings.strength.iter() {
        fs::create_dir_all(format!("Task 1/alpha = {}",alpha))?;

        let v = |r: f64| -> f64 {
            let a = 1.0;
            return if r <= a {
                -alpha.powi(2) / (2.0 * m * a.powi(2))
            } else { 0.0 }
        };

        // Potential output \\
        let mut file = File::create(format!("Task 1/alpha = {}/potential.txt",alpha))?;
        writeln!(file, "{:<3}{:>17}{:>20}", "#", "R [au]", "V(R) [au]")?;
        for i in 0..r.n {
            write!(file, "{:20.12E}", r.points[i])?;
            write!(file, "{:20.12E}", v(r.points[i]))?;
            writeln!(file)?;
        }
        console(&format!("Data successfully written to Task 1/alpha = {}/potential.txt",alpha));

        // Phase-shifts output \\
        let mut file = File::create(format!("Task 1/alpha = {}/phase_shifts.txt",alpha))?;
        let phase_shifts = calculate_phase_shift(&r,v,&e,&ls,m);
        write!(file, "{:<3}{:>17}", "#", "E [au]")?;
        for l in ls.iter() {
            write!(file,"{:>20}",&format!("delta_{} [1]",l))?;
        }
        writeln!(file);
        for i in 0..e.n {
            write!(file, "{:20.12E}", e.points[i])?;
            for j in 0..ls.len() {
                write!(file, "{:20.12E}", phase_shifts[(i,j)])?;
            }
            writeln!(file)?;
        }
        console(&format!("Data successfully written to Task 1/alpha = {}/phase_shifts.txt",alpha));

        // Cross sections output \\
        let mut file = File::create(format!("Task 1/alpha = {}/cross_sections.txt",alpha))?;
        let (cross_sections,cross_section_tot) = calculate_cross_sections(&phase_shifts,&e,&ls);
        write!(file, "{:<3}{:>17}", "#", "E [au]")?;
        for l in ls.iter() {
            write!(file,"{:>20}",&format!("sigma_{} [au^2]",l))?;
        }
        writeln!(file,"{:>20}","sigma_tot [au^2]")?;
        for i in 0..e.n {
            write!(file, "{:20.12E}", e.points[i])?;
            for j in 0..ls.len() {
                write!(file, "{:20.12E}", cross_sections[(i,j)])?;
            }
            write!(file, "{:20.12E}", cross_section_tot[i])?;
            writeln!(file)?;
        }
        console(&format!("Data successfully written to Task 1/alpha = {}/cross_sections.txt",alpha));


        // Limit for zero energy \\
        let phase_shifts = calculate_phase_shift(&r,v,&e_zero,&Vec::from([0]),m);
        let mut file = File::create(format!("Task 1/alpha = {}/zero_energy_limit.txt",alpha))?;
        let (cross_sections,_) = calculate_cross_sections(&phase_shifts,&e_zero,&Vec::from([0]));
        write!(file, "{:<3}{:>17}", "#", "E [au]")?;
        writeln!(file,"{:>20}","sigma_0 [au^2]")?;
        for i in 0..e_zero.n {
            write!(file, "{:20.12E}", e_zero.points[i])?;
            writeln!(file, "{:20.12E}", cross_sections[(i,0)])?;
        }
        console(&format!("Data successfully written to Task 1/alpha = {}/zero_energy_limit.txt",alpha));


        // Bound states calculation \\
        fs::create_dir_all(format!("Task 1/alpha = {}/Matching functions",alpha))?;
        fs::create_dir_all(format!("Task 1/alpha = {}/Eigenfunctions",alpha))?;

        let mut eigenenergies: Vec<DVector<f64>> = Vec::new();
        let mut max_energies = 0;

        for &l in ls.iter() {
            let (energies, functions, calculated_e, matching_f) = calculate_bound_states(&r_bound, v, l, m, settings.bound_ne);
            eigenenergies.push(energies.clone());
            max_energies = max_energies.max(energies.len());

            let mut file = File::create(format!("Task 1/alpha = {}/Matching functions/matching_function_l={}.txt",alpha,l))?;
            writeln!(file, "{:<3}{:>17}{:>20}", "#", "e [au]", "f [a.u.]")?;
            for i in 0..calculated_e.len() {
                write!(file, "{:20.12E}", calculated_e[i])?;
                write!(file, "{:20.12E}", matching_f[i])?;
                writeln!(file)?;
            }
            console(&format!("Data successfully written to Task 1/alpha = {}/Matching functions/matching_function_l={}.txt",alpha,l));

            let mut file = File::create(format!("Task 1/alpha = {}/Eigenfunctions/eigenfunctions_l={}.txt",alpha,l))?;
            write!(file, "{:<3}{:>17}", "#", "R [au]")?;
            for n in 0..energies.len() {
                write!(file,"{:>20}",format!("psi_{} [a.u.]",n))?;
            }
            writeln!(file)?;
            for i in 0..r_bound.n {
                write!(file, "{:20.12E}", r_bound.points[i])?;
                for j in 0..energies.len() {
                    write!(file, "{:20.12E}", functions[(j, i)])?;
                }
                writeln!(file)?;
            }
            console(&format!("Data successfully written to Task 1/alpha = {}/Eigenfunctions/eigenfunctions_l={}.txt",alpha,l));
        }

        let mut file = File::create(format!("Task 1/alpha = {}/bound_state_energies.txt",alpha))?;
        write!(file, "{:<3}{:>1}", "#", "n")?;
        for l in ls.iter() {
            write!(file, "{:>20}",&format!("E_{}", l))?;
        }
        writeln!(file)?;
        for i in 0..max_energies {
            write!(file, "{:4.0}", i)?;
            for energies in eigenenergies.iter() {
                if i < energies.len() {
                    write!(file,"{:20.12E}",energies[i])?;
                }
                else {
                    write!(file,"{:>20}","NaN")?;
                }
            }
            writeln!(file)?;
        }
        console(&format!("Data successfully written to Task 1/alpha = {}/bound_state_energies.txt",alpha));
    }






    // Task 2 Output \\
    let r = Grid::new_n(settings.rmin2,settings.rmax2,settings.nr2);
    let r_bound = Grid::new_n(settings.rmin2,settings.rmax2 * settings.bound_mult as f64,settings.nr2 * settings.bound_mult);
    let e = Grid::new_d(settings.emin2,settings.emax2,settings.de2);
    let e_zero = Grid::new_log(std::f64::EPSILON,settings.de2,settings.zero_limit);
    let ls: Vec<usize> = settings.ls2;
    let m = settings.mass2;

    let v = |r: f64| -> f64 {
        -3.0 * (-r.powi(2)/4.0).exp() + (-(r-3.0).powi(2)).exp()
    };

    // Potential output \\
    let mut file = File::create("Task 2/potential.txt")?;
    writeln!(file, "{:<3}{:>17}{:>20}", "#", "R [au]", "V(R) [au]")?;
    for i in 0..r.n {
        write!(file, "{:20.12E}", r.points[i])?;
        write!(file, "{:20.12E}", v(r.points[i]))?;
        writeln!(file)?;
    }
    console("Data successfully written to Task 2/potential.txt");

    // Phase-shifts output \\
    let mut file = File::create("Task 2/phase_shifts.txt")?;
    let phase_shifts = calculate_phase_shift(&r,v,&e,&ls,m);
    write!(file, "{:<3}{:>17}", "#", "E [au]")?;
    for l in ls.iter() {
        write!(file,"{:>20}",&format!("delta_{} [1]",l))?;
    }
    writeln!(file);
    for i in 0..e.n {
        write!(file, "{:20.12E}", e.points[i])?;
        for j in 0..ls.len() {
            write!(file, "{:20.12E}", phase_shifts[(i,j)])?;
        }
        writeln!(file)?;
    }
    console("Data successfully written to Task 2/phase_shifts.txt");

    // Cross sections output \\
    let mut file = File::create("Task 2/cross_sections.txt")?;
    let (cross_sections,cross_section_tot) = calculate_cross_sections(&phase_shifts,&e,&ls);
    write!(file, "{:<3}{:>17}", "#", "E [au]")?;
    for l in ls.iter() {
        write!(file,"{:>20}",&format!("sigma_{} [au^2]",l))?;
    }
    writeln!(file,"{:>20}","sigma_tot [au^2]")?;
    for i in 0..e.n {
        write!(file, "{:20.12E}", e.points[i])?;
        for j in 0..ls.len() {
            write!(file, "{:20.12E}", cross_sections[(i,j)])?;
        }
        write!(file, "{:20.12E}", cross_section_tot[i])?;
        writeln!(file)?;
    }
    console("Data successfully written to Task 2/cross_sections.txt");


    // Limit for zero energy \\
    let phase_shifts = calculate_phase_shift(&r,v,&e_zero,&Vec::from([0]),m);
    let mut file = File::create("Task 2/zero_energy_limit.txt")?;
    let (cross_sections,_) = calculate_cross_sections(&phase_shifts,&e_zero,&Vec::from([0]));
    write!(file, "{:<3}{:>17}", "#", "E [au]")?;
    writeln!(file,"{:>20}","sigma_0 [au^2]")?;
    for i in 0..e_zero.n {
        write!(file, "{:20.12E}", e_zero.points[i])?;
        writeln!(file, "{:20.12E}", cross_sections[(i,0)])?;
    }
    console("Data successfully written to Task 2/zero_energy_limit.txt");


    // Bound states calculation \\
    fs::create_dir_all("Task 2/Matching functions")?;
    fs::create_dir_all("Task 2/Eigenfunctions")?;

    let mut eigenenergies: Vec<DVector<f64>> = Vec::new();
    let mut max_energies = 0;

    for &l in ls.iter() {
        let (energies, functions, calculated_e, matching_f) = calculate_bound_states(&r_bound, v, l, m, settings.bound_ne);
        eigenenergies.push(energies.clone());
        max_energies = max_energies.max(energies.len());

        let mut file = File::create(format!("Task 2/Matching functions/matching_function_l={}.txt",l))?;
        writeln!(file, "{:<3}{:>17}{:>20}", "#", "e [au]", "f [a.u.]")?;
        for i in 0..calculated_e.len() {
            write!(file, "{:20.12E}", calculated_e[i])?;
            write!(file, "{:20.12E}", matching_f[i])?;
            writeln!(file)?;
        }
        console(&format!("Data successfully written to Task 2/Matching functions/matching_function_l={}.txt",l));

        let mut file = File::create(format!("Task 2/Eigenfunctions/eigenfunctions_l={}.txt",l))?;
        write!(file, "{:<3}{:>17}", "#", "R [au]")?;
        for n in 0..energies.len() {
            write!(file,"{:>20}",format!("psi_{} [a.u.]",n))?;
        }
        writeln!(file)?;
        for i in 0..r_bound.n {
            write!(file, "{:20.12E}", r_bound.points[i])?;
            for j in 0..energies.len() {
                write!(file, "{:20.12E}", functions[(j, i)])?;
            }
            writeln!(file)?;
        }
        console(&format!("Data successfully written to Task 2/Eigenfunctions/eigenfunctions_l={}.txt",l));
    }

    let mut file = File::create("Task 2/bound_state_energies.txt")?;
    write!(file, "{:<3}{:>1}", "#", "n")?;
    for l in ls.iter() {
        write!(file, "{:>20}",&format!("E_{}", l))?;
    }
    writeln!(file)?;
    for i in 0..max_energies {
        write!(file, "{:4.0}", i)?;
        for energies in eigenenergies.iter() {
            if i < energies.len() {
                write!(file,"{:20.12E}",energies[i])?;
            }
            else {
                write!(file,"{:>20}","NaN")?;
            }
        }
        writeln!(file)?;
    }
    console("Data successfully written to Task 2/bound_state_energies.txt");

    Ok(())
}
