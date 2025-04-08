use core::panic;
use std::ops::{Add, Mul, Sub};

use eframe::egui::{Color32, Vec2};
use rand::Rng;
use rand_distr::{Distribution, Normal};

use super::NanError;

/// Linear interpolation between two variables.
///
/// When `r` is 0, the output will equal to `f1`. When `r` is 1, the output will equal to `f2`.
///
/// And when `r` is in between, the output linearly varies from `f1` to `f2` as `r` increases.
///
/// The type of `f1` and `f2` can be any type that implements all [Add], [Sub], and [Mul<f32>].
///
/// # Examples
///
/// ```
/// use thermal_model::linear_interp;
/// assert_eq!(linear_interp(1.0, 2.0, 0.5), 1.5);
/// assert_eq!(linear_interp(1.0, 2.0, 0.1), 1.1);
/// assert_eq!(linear_interp(1.0, 2.0, 0.9), 1.9);
/// ```
pub fn linear_interp<F: Add<Output = F> + Sub<Output = F> + Mul<f32, Output = F> + Clone>(
    f1: F,
    f2: F,
    r: f32,
) -> F {
    f1.clone() + (f2 - f1) * r
}

pub fn color_interp(c1: Color32, c2: Color32, ratio: f32) -> Color32 {
    let r1 = c1.r() as f32;
    let g1 = c1.g() as f32;
    let b1 = c1.b() as f32;
    let r2 = c2.r() as f32;
    let g2 = c2.g() as f32;
    let b2 = c2.b() as f32;
    let r = (r1 + (r2 - r1) * ratio) as u8;
    let g = (g1 + (g2 - g1) * ratio) as u8;
    let b = (b1 + (b2 - b1) * ratio) as u8;
    Color32::from_rgb(r, g, b)
}

pub fn sample_rand_velocity(rng: &mut impl Rng, stdev: f32) -> Vec2 {
    let normal = Normal::new(0.0, stdev).unwrap();
    Vec2::new(normal.sample(rng), normal.sample(rng))
}

/// Cross product of two 2D vectors.
/// 
/// # Examples
/// 
/// ```
/// use eframe::egui::Vec2;
/// use thermal_model::cross;
/// let v1 = Vec2::new(1.0, 2.0);
/// let v2 = Vec2::new(3.0, 4.0);
/// assert_eq!(cross(v1, v2), -2.0);
/// ```
pub fn cross(v1: Vec2, v2: Vec2) -> f32 {
    v1.x * v2.y - v1.y * v2.x
}

/// Solve the equation set
/// 
/// $ax+by=c$
/// 
/// $ax^2+by^2=d$
/// 
/// The function returns the solution $(x, y)$.
pub fn solve_conservation_equation(a: f64, b: f64, c: f64, d: f64, x_positive: bool) -> Result<(f64, f64), NanError> {
    let delta = a * b * (d * (a + b) - c * c);
    if delta < 0.0 {
        // return Err(NanError::NegativeSquareRoot);
        return Ok((c / (a + b), c / (a + b)));
    }
    if x_positive {
        Ok((
            (a * c + delta.sqrt()) / (a * (a + b)),
            (b * c - delta.sqrt()) / (b * (a + b)),
        ))
    } else {
        Ok((
            (a * c - delta.sqrt()) / (a * (a + b)),
            (b * c + delta.sqrt()) / (b * (a + b)),
        ))
    }
}