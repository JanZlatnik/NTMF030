/*

MATH MODULE

Contains: Basic mathematical functions and operations

Last revision: 26/11/2024 by JZ

*/
use std::fmt::Debug;
use nalgebra::{Complex, DMatrix, DVector};
use num::traits::{Num, Zero};
use std::ops::{Add, Mul};
use num_traits::{Float};
use std::collections::HashMap;
use num::complex::ComplexFloat;

// Grid definition \\
#[derive(Clone)]
pub struct Grid {
    pub points: Vec<f64>,
    pub d:      f64,
    pub n:      usize,
    pub start:  f64,
    pub end:    f64,
}

impl Grid {
    pub(crate) fn new_n(start: f64, end: f64, n: usize) -> Grid {
        let d = (end - start)/((n - 1) as f64);
        let points = (0.. n)
            .map(|i| start + i as f64 * d)
            .collect();

        Grid{
            points,
            d,
            n,
            start,
            end
        }
    }

    pub(crate) fn new_d(start: f64, end: f64, spacing: f64) -> Grid {
        let d = spacing;
        let n = ((end - start) / d).ceil() as usize + 1;

        let points = (0..n)
            .map(|i| start + i as f64 * d)
            .collect();


        Grid {
            points,
            d,
            n,
            start,
            end,
        }
    }

    pub(crate) fn tail(&mut self, t:f64) {
        if t <= 1.0 {
            return;
        }
        let new_n = ((self.n as f64) * t).ceil() as usize;
        for i in 0..(new_n-self.n) {
            self.points.push(self.points[self.n+i] + self.d);
        }
        self.end = self.points[new_n-1];
        self.n = self.points.len();
    }
}




// Integration \\
pub trait Integrable<T> {
    fn integrate(&self, dx: f64) -> T;
    fn square_integrate(&self, dx: f64) -> T;
    fn primitive(&self, dx: f64) -> DVector<T>;
}

impl<T: Num + Copy + Add<Output=T> + Mul<f64, Output=T> + Zero + Default + Debug + 'static> Integrable<T> for DVector<T> {
    fn integrate(&self, dx: f64) -> T {
        let n = self.len();
        if n < 2 {
            panic!("[Integrate Error]: Not enough elements")
        }
        let mut sum = T::zero();
        for i in 1..n-1 {
            sum = sum + self[i];
        }
        sum = sum + (self[0]+self[n-1]) * 0.5;
        sum * dx
    }

    fn square_integrate(&self, dx: f64) -> T {
        let n = self.len();
        if n < 2 {
            panic!("[Integrate Error]: Not enough elements")
        }
        let mut sum = T::zero();
        for i in 1..n-1 {
            sum = sum + self[i] * self[i];
        }
        sum = sum + (self[0] * self[0] +self[n-1] * self[n-1]) * 0.5;
        sum * dx
    }


    fn primitive(&self, dx: f64) -> DVector<T> {
        let n = self.len();
        let mut primitive_vec =  DVector::from_vec(vec![T::default(); n]);
        if n < 2 {
            panic!("[Primitive Error]: Not enough elements")
        }
        for i in 1..n {
            primitive_vec[i] = primitive_vec[i-1] + (self[i-1]+self[i])* 0.5 * dx;
        }
        primitive_vec
    }
}



// Riccati-Bessel and Riccati-Neumann function definitions for real arguments \\
pub fn riccati_bessel(l: usize, x: f64) -> f64 {
    match l {
        0 => x.sin(),
        1 => x.sin() / x - x.cos(),
        _ => {
            let mut j0 = x.sin();
            let mut j1 = x.sin() / x - x.cos();
            let mut jl = 0.0;
            for i in 2..=l {
                jl = ((2 * i - 1) as f64 / x) * j1 - j0;
                j0 = j1;
                j1 = jl;
            }
            jl
        }
    }
}

pub fn riccati_neumann(l: usize, x: f64) -> f64 {
    match l {
        0 => x.cos(),
        1 => x.cos() / x + x.sin(),
        _ => {
            let mut n0 = x.cos();
            let mut n1 = x.cos() / x + x.sin();
            let mut nl = 0.0;
            for i in 2..=l {
                nl = ((2 * i - 1) as f64 / x) * n1 - n0;
                n0 = n1;
                n1 = nl;
            }
            nl
        }
    }
}