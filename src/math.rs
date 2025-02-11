/*

MATH MODULE

Contains: Basic mathematical functions and operations

Last revision: 07/02/2025 by JZ

*/
use std::fmt::Debug;
use nalgebra::{Complex, DVector};
use num::traits::{Num, Zero};
use std::ops::{Add, Mul};
use num_traits::{Float};
use crate::core::console;

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

    pub(crate) fn new_log(start: f64, end: f64, n: usize) -> Grid {
        let log_start = start.ln();
        let log_end = end.ln();
        let d = (log_end - log_start) / ((n - 1) as f64);

        let points = (0..n)
            .map(|i| (log_start + i as f64 * d).exp())
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
        if t == 1.0 {
            return;
        }
        else if t < 0.0 {
            console("[ERROR]: tail in grid.tail cannot be smaller than 0");
            return;
        }
        else if t >  1.0 {
            let new_n = ((self.n as f64) * t).ceil() as usize;
            for i in 0..(new_n-self.n) {
                self.points.push(self.points[self.n-1+i] + self.d);
            }
            self.end = self.points[new_n-1];
            self.n = self.points.len();
        }
        else if t < 1.0 {
            let new_n = ((self.n as f64) * t).ceil() as usize;
            self.points.truncate(new_n);
            self.end = self.points[new_n-1];
            self.n = self.points.len();
        }
    }
}




// Integration \\
#[derive(Clone, Copy)]
pub enum INTMethod {
    Trapezoid,
    Simpson,
    INT4
}
#[derive(Clone, Copy)]
pub enum PrimitiveMethod {
    P2,
    P4
}
pub trait Integrable<T> {
    fn integrate(&self, dx: f64, method: INTMethod) -> T;
    fn square_integrate(&self, dx: f64, method: INTMethod) -> f64;
    fn primitive(&self, dx: f64, method: PrimitiveMethod) -> DVector<T>;
}

trait ComplexOrReal {
    fn abs(&self) -> f64;
}
impl ComplexOrReal for f64 {
    fn abs(&self) -> f64 { num_traits::Float::abs(*self)  }
}
impl<F: Float> ComplexOrReal for Complex<F> {
    fn abs(&self) -> f64 { self.norm().to_f64().unwrap_or(0.0) }
}

impl<T: Num + Copy + Add<Output=T> + Mul<f64, Output=T> + Zero + Default + Debug + 'static> Integrable<T> for DVector<T>
where T: ComplexOrReal {
    fn integrate(&self, dx: f64, method: INTMethod) -> T {
        let n = self.len();
        let mut sum = T::zero();
        match method {
            INTMethod::Trapezoid => {
                if n < 2 {
                    panic!("[Integrate Error]: Not enough elements")
                }
                for i in 1..n-1 {
                    sum = sum + self[i];
                }
                sum = sum + (self[0]+self[n-1]) * 0.5;
            }
            INTMethod::Simpson => {
                if n < 3 {
                    panic!("[Integrate Error]: Not enough elements")
                }
                let m = n - 1 + (n % 2);
                sum = sum + self[0] + self[m-1];
                for i in 1..=(m - 1) / 2 {
                    sum = sum + self[2*i-1] * 4.0
                }
                for i in 1..=(m-1)/2 - 1 {
                    sum = sum + self[2*i] * 2.0
                }
                sum = sum * (1.0 / 3.0);
                if m != n {
                    sum = sum + self[n-1] * (3.0/8.0) + self[n-2] * (19.0/24.0) - self[n-3] * (5.0/24.0) + self[n-4] * (1.0/ 24.0);
                }
            }
            INTMethod::INT4 => {
                if n < 6 {
                    panic!("[Integrate Error]: Not enough elements")
                }
                for i in 3..n-3 {
                    sum = sum + self[i]
                }
                sum = sum + (self[0]+self[n-1]) * (3.0/8.0) + (self[1]+self[n-2]) * (7.0 / 6.0) + (self[2]+self[n-3]) * (23.0 / 24.0)
            }
        }

        sum * dx
    }

    fn square_integrate(&self, dx: f64, method: INTMethod) -> f64 {
        let squared_vector = DVector::from_vec(
            self.iter().map(|&x| x.abs() * x.abs()).collect()
        );
        squared_vector.integrate(dx, method)
    }


    fn primitive(&self, dx: f64, method: PrimitiveMethod) -> DVector<T> {
        let n = self.len();
        let mut primitive_vec =  DVector::from_vec(vec![T::zero(); n]);
        match method {
            PrimitiveMethod::P2 => {
                if n < 2 {
                    panic!("[Primitive Error]: Not enough elements")
                }
                for i in 1..n {
                    primitive_vec[i] = primitive_vec[i-1] + (self[i-1]+self[i])* 0.5  * dx;
                }
            }
            PrimitiveMethod::P4 => {
                if n < 4 {
                    panic!("[Primitive Error]: Not enough elements")
                }
                primitive_vec[1] = (self[0] * (3.0/8.0)  + self[1] * (19.0/24.0) - self[2] * (5.0/24.0) + self[3] * (1.0/24.0))  * dx;
                for i in 2..n-1 {
                    primitive_vec[i] = primitive_vec[i-1] + ((self[i-1]+self[i]) * (13.0/24.0) - (self[i-2]+self[i+1]) * (1.0/24.0)) * dx;
                }
                primitive_vec[n-1] = primitive_vec[n-2] + (self[n-1] * (3.0/8.0) + self[n-2] * (19.0/24.0) - self[n-3] * (5.0/24.0) + self[n-4] * (1.0/24.0)) * dx;
            }

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