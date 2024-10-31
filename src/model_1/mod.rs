use eframe::egui::Vec2;
use rand::Rng;

use crate::{
    model::{Model, Molecule},
    utils, MoleculeType,
};

/// Naive implementation of gas model
///
/// Use O(N^2) algorithm to check and process every collision within every frame.
pub struct Model1<T: MoleculeType> {
    molecules: Vec<Molecule<T>>,
    width: f32,
    height: f32,
}

impl<T: MoleculeType + Default> Model1<T> {
    pub fn new(width: f32, height: f32, radius: f32, mass: f32, num_molecule: usize) -> Self {
        let mol_type = T::default();
        let mut rng = rand::thread_rng();
        let molecules = (0..num_molecule)
            .map(|_| Molecule {
                pos: Vec2::new(
                    rng.gen_range(radius..(width - radius)),
                    rng.gen_range(radius..(height - radius)),
                ),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                mol_type,
            })
            .collect();
        Self {
            molecules,
            width,
            height,
        }
    }
}

impl<T: MoleculeType> Model1<T> {
    /// Handle collision between molecules.
    ///
    /// Returns the total number of collisions occured, including both molecule-molecule and molecule-wall collisions.
    fn handle_collision(&mut self) -> usize {
        let mut total_collisions = 0;
        // handle collisions between molecules
        let num_molecules = self.molecules.len();
        for i in 0..num_molecules {
            let (left, right) = self.molecules.split_at_mut(i + 1);
            let mol1 = left[i];
            for j in 0..(num_molecules - i - 1) {
                let mol2 = right[j];
                if Vec2::length(mol1.pos - mol2.pos)
                    <= mol1.mol_type.radius() + mol2.mol_type.radius()
                {
                    utils::collision_2_molecules(&mut left[i], &mut right[j]);
                    total_collisions += 1;
                }
            }
        }
        // handle collisions with walls
        for m in self.molecules.iter_mut() {
            let radius = m.mol_type.radius();
            if m.pos.x - radius < 0.0 {
                m.pos.x = 2.0 * radius - m.pos.x;
                m.vel.x = f32::abs(m.vel.x);
                total_collisions += 1;
            }
            if m.pos.x + radius > self.width {
                m.pos.x = 2.0 * self.width - 2.0 * radius - m.pos.x;
                m.vel.x = -f32::abs(m.vel.x);
                total_collisions += 1;
            }
            if m.pos.y - radius < 0.0 {
                m.pos.y = 2.0 * radius - m.pos.y;
                m.vel.y = f32::abs(m.vel.y);
                total_collisions += 1;
            }
            if m.pos.y + radius > self.height {
                m.pos.y = 2.0 * self.height - 2.0 * radius - m.pos.y;
                m.vel.y = -f32::abs(m.vel.y);
                total_collisions += 1;
            }
        }
        total_collisions
    }
}

impl<T: MoleculeType> Model for Model1<T> {
    type Type = T;
    type AdvanceReturnType = usize;

    fn construct(
        width: f32,
        height: f32,
        num_molecule: usize,
        constructor: impl FnMut(usize) -> Molecule<Self::Type>,
    ) -> Self {
        let molecules = (0..num_molecule).map(constructor).collect();
        Self {
            molecules,
            width,
            height,
        }
    }

    fn dimension(&self) -> (f32, f32) {
        (self.width, self.height)
    }

    fn num_molecule(&self) -> usize {
        self.molecules.len()
    }

    fn get_molecules(&self) -> impl Iterator<Item = Molecule<T>> {
        self.molecules.iter().cloned()
    }

    fn advance(&mut self, dt: f32) -> usize {
        // move all molecules along their speeds
        for m in self.molecules.iter_mut() {
            m.pos += m.vel * dt;
        }
        // handle collisions multiple times
        for _ in 0..4 {
            self.handle_collision();
        }
        self.handle_collision() // Return the last collision count. Expected to be 0.
    }
}
