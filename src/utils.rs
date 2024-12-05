use std::ops::{Add, Mul, Sub};

use eframe::egui::{Color32, Vec2};
use rand::Rng;
use rand_distr::{Distribution, Normal};

use crate::model::{Molecule, MoleculeType};

/// The impulse acted on the first molecule during the collision
///
/// The impulse on the second molecule is the negative of this value.
#[inline]
pub fn normal_impulse(m1: f32, m2: f32, v1: f32, v2: f32) -> f32 {
    2.0 * m1 * m2 * (v2 - v1) / (m1 + m2)
}

const TANGENTIAL_IMPULSE_CONSTANT: f32 = 0.2;

pub fn tangential_impulse(m1: f32, m2: f32, v1: f32, v2: f32, w1: f32, w2: f32) -> f32 {
    // assume that the impulse on tangential direction is proportional to the difference in angular
    // velocity
    -TANGENTIAL_IMPULSE_CONSTANT * (w1 + w2) // since the two molecules spins at opposite direction at the contact point
}

pub fn collision_2_molecules<T: MoleculeType>(mol1: &mut Molecule<T>, mol2: &mut Molecule<T>) {
    // letters
    // m1, m2: mass
    // v1, v2: velocity
    // r1, r2: radius
    // w1, w2: angular velocity
    // d_v1, d_v2: velocity change
    // d_v1_n, d_v2_n: velocity change on normal direction
    // d_v1_t, d_v2_t: velocity change on tangential direction

    // calculate velocities after collision
    let normal = (mol2.pos - mol1.pos).normalized();
    let m1 = mol1.mol_type.mass();
    let m2 = mol2.mol_type.mass();
    let v1 = mol1.vel.dot(normal);
    let v2 = mol2.vel.dot(normal);
    let r1 = mol1.mol_type.radius();
    let r2 = mol2.mol_type.radius();
    let impulse_n = normal_impulse(m1, m2, v1, v2);
    let mut d_v1 = impulse_n / m1 * normal;
    let mut d_v2 = -impulse_n / m2 * normal;

    match &mut (mol1.orient, mol2.orient) {
        (Some((_, w1)), Some((_, w2))) => {
            // calculate change in velocity on tangential direction
            let tangential = -normal.rot90();
            let impulse_t = tangential_impulse(m1, m2, v1, v2, *w1, *w2);
            d_v1 += impulse_t / m1 * tangential;
            d_v2 -= impulse_t / m2 * tangential;

            // change angular velocity
            *w1 += r1 * impulse_t / mol1.mol_type.inertia();
            *w2 -= r2 * impulse_t / mol2.mol_type.inertia();
        }
        _ => {}
    }

    // update the two molecules' velocities
    mol1.vel += d_v1;
    mol2.vel += d_v2;

    // update the two molecules' positions to avoid overlap
    let update_dist = (r1 + r2 - Vec2::length(mol1.pos - mol2.pos)) / 2.0;
    mol1.pos -= normal * update_dist;
    mol2.pos += normal * update_dist;
}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collision() {
        assert_eq!(normal_impulse(1.0, 1.0, 1.0, 0.0), -1.0);
        assert_eq!(normal_impulse(1.0, 1.0, 1.0, -1.0), -2.0);
    }

    #[test]
    fn test_interpolation() {
        assert_eq!(linear_interp(1.0, 2.0, 0.5), 1.5);
        assert_eq!(linear_interp(1.0, 2.0, 0.1), 1.1);
        assert_eq!(linear_interp(1.0, 2.0, 0.9), 1.9);
    }
}
