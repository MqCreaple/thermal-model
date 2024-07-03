use eframe::egui::Vec2;
use rand::Rng;

use crate::{model::{Model, Molecule}, utils};

pub mod kd_tree;

/// Naive implementation of gas model
/// 
/// Use O(N^2) algorithm to check and process every collision within every frame.
pub struct Model2 {
    molecules: Vec<Molecule>,
    width: f32,
    height: f32,
    kd_tree_max_elem: usize,
}

#[derive(Clone, Copy, Debug)]
struct IndexedMolecule {
    index: usize,
    molecule: Molecule,
}

impl kd_tree::PositionKey for IndexedMolecule {
    fn get_position(&self) -> Vec2 {
        self.molecule.pos
    }
}

impl Model2 {
    pub fn new(width: f32, height: f32, radius: f32, mass: f32, num_molecule: usize, kd_tree_max_elem: usize) -> Self {
        let mut rng = rand::thread_rng();
        let molecules = (0..num_molecule).map(|_| {
            Molecule {
                pos: Vec2::new(rng.gen_range(radius..(width - radius)), rng.gen_range(radius..(height - radius))),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                radius, mass
            }
        }).collect();
        Self {
            molecules, width, height, kd_tree_max_elem
        }
    }

    /// Handle collision between molecules.
    /// 
    /// Returns the total number of collisions occured, including both molecule-molecule and molecule-wall collisions.
    fn handle_collision(&mut self) -> usize {
        let mut total_collisions = 0;
        // build a KD tree of all molecules
        let mut indexed_molecule = self.molecules.iter().enumerate().map(|(i, m)| IndexedMolecule { index: i, molecule: *m }).collect::<Vec<_>>();
        let kd_tree = kd_tree::Node::new(&mut indexed_molecule, self.kd_tree_max_elem, kd_tree::KDNodeDivideBy::X);
        // handle collisions between molecules
        let num_molecules = self.molecules.len();
        for i in 0..num_molecules {
            let mol1 = self.molecules[i];
            let neighbors = kd_tree.query_circle(mol1.pos, mol1.radius * 2.0);
            for indexed_mol2 in neighbors {
                let j = indexed_mol2.index;
                let mol2 = &self.molecules[indexed_mol2.index];
                if i != j && Vec2::length(mol1.pos - mol2.pos) < mol1.radius + mol2.radius {
                    let (slice_1, slice_2) = self.molecules.split_at_mut(i.max(j));
                    // handle collision
                    utils::collision_2_molecules(&mut slice_1[i.min(j)], &mut slice_2[0]);
                    total_collisions += 1;
                }
                // TODO: bug: when the position of the molecule is updated here, the KD tree is not updated
                // Therefore, the collisions calculated with the KD tree might not be accurate.
            }
        }
        // handle collisions with walls
        for m in self.molecules.iter_mut() {
            if m.pos.x - m.radius < 0.0 {
                m.pos.x = 2.0 * m.radius - m.pos.x;
                m.vel.x = f32::abs(m.vel.x);
                total_collisions += 1;
            }
            if m.pos.x + m.radius > self.width {
                m.pos.x = 2.0 * self.width - 2.0 * m.radius - m.pos.x;
                m.vel.x = -f32::abs(m.vel.x);
                total_collisions += 1;
            }
            if m.pos.y - m.radius < 0.0 {
                m.pos.y = 2.0 * m.radius - m.pos.y;
                m.vel.y = f32::abs(m.vel.y);
                total_collisions += 1;
            }
            if m.pos.y + m.radius > self.height {
                m.pos.y = 2.0 * self.height - 2.0 * m.radius - m.pos.y;
                m.vel.y = -f32::abs(m.vel.y);
                total_collisions += 1;
            }
        }
        total_collisions
    }
}

impl Model for Model2 {
    type AdvanceReturnType = usize;

    fn dimension(&self) -> (f32, f32) {
        (self.width, self.height)
    }

    fn num_molecule(&self) -> usize {
        self.molecules.len()
    }

    fn get_molecules(&self) -> impl Iterator<Item = Molecule> {
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
        self.handle_collision()
    }
}