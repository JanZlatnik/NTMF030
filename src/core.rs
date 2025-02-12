/*

CORE MODULE

Contains: Basic mathematical and physical constants, settings initialization, grid definition

Last revision: 26/11/2024 by JZ

*/
use std::path::Path;
use chrono::Local;
use serde::Deserialize;

// Useful physical & mathematical constants \\
pub const PI:       f64 = std::f64::consts::PI;         // Pi
pub const EULER:    f64 = std::f64::consts::E;          // Euler's number
pub const PHYS_A0:  f64 = 0.52917706;                   // Bohr radius [angstroms]
pub const PHYS_H0:  f64 = 27.2113961;                   // 1 Hartree in [eV]
pub const PHYS_K:   f64 = 8.617385e-5;                  // Boltzmann constant in [eV/K]
pub const PHYS_ME:  f64 = 9.1093897e-31;                // electron mass [kg]
pub const PHYS_U:   f64 = 1.6605402e-27;                // (unified) atomic mass unit [kg]
pub const HBAR:     f64 = 1.0;                          // Planck's constant [au]



// Console message with time print \\
pub fn console(message: &str) {
    let now = Local::now();
    let time = now.format("%H:%M:%S").to_string();
    println!("[{}]: {}", time, message);
}



// Settings initialization \\
#[derive(Deserialize,Debug)]
pub struct Settings {
    // Grid settings Task 1
    pub rmin1: f64,
    pub rmax1: f64,
    pub nr1: usize,

    // Energy settings Task 1
    pub emin1: f64,
    pub emax1: f64,
    pub de1: f64,

    // Grid settings Task 2
    pub rmin2: f64,
    pub rmax2: f64,
    pub nr2: usize,

    // Energy settings Task 2
    pub emin2: f64,
    pub emax2: f64,
    pub de2: f64,

    // Calculation settings Task 1
    pub ls1: Vec<usize>,
    pub mass1: f64,
    pub strength: Vec<f64>,

    // Calculation settings Task 2
    pub ls2: Vec<usize>,
    pub mass2: f64,

    // Bound states settings
    pub bound_ne: usize,
    pub bound_mult: usize,

    // Other settings
    pub zero_limit: usize,
}

impl Settings {
    pub(crate) fn initialize<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        // read settings
        let content = std::fs::read_to_string(path)
            .map_err(|e| format!("Failed to read settings file: {}", e))?;
        let settings: Settings = toml::from_str(&content)
            .map_err(|e| format!("Failed to parse TOML: {}", e))?;

        Ok(settings)
    }
}


