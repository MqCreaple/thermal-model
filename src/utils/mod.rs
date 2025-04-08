mod math;
mod error;

use eframe::egui::Vec2;
use rand::Rng;

use crate::model::{Molecule, MoleculeType};

pub use math::*;
pub use error::*;

/// The impulse acted on the first molecule during the collision
///
/// The impulse on the second molecule is the negative of this value.
#[inline]
pub fn normal_impulse(m1: f32, m2: f32, v1_normal: f32, v2_normal: f32) -> f32 {
    2.0 * m1 * m2 * (v2_normal - v1_normal) / (m1 + m2)
}

const TANGENTIAL_IMPULSE_CONSTANT: f32 = 0.05;

#[inline]
pub fn tangential_impulse(m1: f32, m2: f32, v1: f32, v2: f32, w1: f32, w2: f32) -> f32 {
    // assume that the impulse on tangential direction is proportional to the difference in angular
    // velocity
    -TANGENTIAL_IMPULSE_CONSTANT * (w1 + w2) // since the two molecules spins at opposite direction at the contact point
}

const LEARN_RATE: f32 = 0.01;

/// Use gradient descent to correct the velocity of mol1 and mol2 to make their total translational
/// kinetic energy equal to `energy`.
#[inline]
pub fn velocity_correction<T: MoleculeType>(mol1: &mut Molecule<T>, mol2: &mut Molecule<T>, energy: f32) {
    let m1 = mol1.mol_type.mass();
    let m2 = mol2.mol_type.mass();
    for _ in 0..3 {
        let cur_energy = 0.5 * (m1 * mol1.vel.length_sq() + m2 * mol2.vel.length_sq());
        let momentum = LEARN_RATE * 2.0 * (energy - cur_energy) * (mol1.vel - mol2.vel);
        mol1.vel += momentum / m1;
        mol2.vel -= momentum / m2;
    }
}

pub fn collision_2_molecules<T: MoleculeType>(mol1: &mut Molecule<T>, mol2: &mut Molecule<T>, rng: &mut impl Rng) {
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
    let energy = 0.5 * (m1 * mol1.vel.length_sq() + m2 * mol2.vel.length_sq());
    let v1_normal = mol1.vel.dot(normal);
    let v2_normal = mol2.vel.dot(normal);
    let r1 = mol1.mol_type.radius();
    let r2 = mol2.mol_type.radius();
    let impulse_n = normal_impulse(m1, m2, v1_normal, v2_normal);
    let mut d_v1 = impulse_n / m1 * normal;
    let mut d_v2 = -impulse_n / m2 * normal;
    
    // update the two molecules' velocities
    mol1.vel += d_v1;
    mol2.vel += d_v2;
    velocity_correction(mol1, mol2, energy);
    // TODO: update the angular velocity

    // update the two molecules' positions to avoid overlap
    let update_dist = (r1 + r2 - Vec2::length(mol1.pos - mol2.pos)) / 2.0;
    mol1.pos -= normal * update_dist;
    mol2.pos += normal * update_dist;
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

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

    #[test]
    fn test_conservation_equation() {
        let (x, y) = solve_conservation_equation(1.0, 1.0, 1.0, 1.0, true).unwrap();
        assert_relative_eq!(x + y, 1.0);
        assert_relative_eq!(x * x + y * y, 1.0);

        let (x, y) = solve_conservation_equation(1.0, 1.0, 1.0, 1.0, false).unwrap();
        assert_relative_eq!(x + y, 1.0);
        assert_relative_eq!(x * x + y * y, 1.0);

        let (x, y) = solve_conservation_equation(1.0, 2.0, 3.0, 4.0, true).unwrap();
        assert_relative_eq!(x + 2.0 * y, 3.0);
        assert_relative_eq!(x * x + 2.0 * y * y, 4.0);

        let (x, y) = solve_conservation_equation(1.0, 2.0, 3.0, 4.0, false).unwrap();
        assert_relative_eq!(x + 2.0 * y, 3.0);
        assert_relative_eq!(x * x + 2.0 * y * y, 4.0);

        let (x, y) = solve_conservation_equation(1.0, -2.0, 3.0, -4.0, true).unwrap();
        assert_relative_eq!(x - 2.0 * y, 3.0);
        assert_relative_eq!(x * x - 2.0 * y * y, -4.0);

        let (x, y) = solve_conservation_equation(1.0, -2.0, 3.0, -4.0, false).unwrap();
        assert_relative_eq!(x - 2.0 * y, 3.0);
        assert_relative_eq!(x * x - 2.0 * y * y, -4.0);
    }
}
