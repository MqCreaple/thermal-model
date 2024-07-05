use std::collections::BTreeSet;

use eframe::egui::Vec2;
use num::Zero;
use ordered_float::NotNan;
use rand::Rng;

use crate::{model::{Model, Molecule, MoleculeType}, utils};

pub mod kd_tree;

/// Naive implementation of gas model
/// 
/// Use O(N^2) algorithm to check and process every collision within every frame.
pub struct Model2<T: MoleculeType> {
    molecules: Vec<Molecule<T>>,
    width: f32,
    height: f32,
    kd_tree_max_elem: usize,
}

/// A molecule with its index in the molecule array
#[derive(Clone, Copy, Debug)]
struct IndexedMolecule {
    index: usize,
    pos: Vec2,
    radius: f32,
}

impl kd_tree::PositionKey for IndexedMolecule {
    fn get_position(&self) -> Vec2 {
        self.pos
    }
}

#[derive(PartialEq, Eq, Debug)]
struct CollisionEntry(usize, usize, NotNan<f32>);

impl PartialOrd for CollisionEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.2.partial_cmp(&other.2)
    }
}

impl Ord for CollisionEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.2.cmp(&other.2)
    }
}

impl<T: MoleculeType + Default> Model2<T> {
    pub fn new(width: f32, height: f32, num_molecule: usize, kd_tree_max_elem: usize) -> Self {
        let default = T::default();
        let radius = default.radius();
        let mut rng = rand::thread_rng();
        let molecules = (0..num_molecule).map(|_| {
            Molecule {
                pos: Vec2::new(rng.gen_range(radius..(width - radius)), rng.gen_range(radius..(height - radius))),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                mol_type: default,
            }
        }).collect();
        Self {
            molecules, width, height, kd_tree_max_elem,
        }
    }
}

impl<T: MoleculeType> Model2<T> {
    /// Handle collision between molecules.
    /// 
    /// Returns the total number of collisions occured, including both molecule-molecule and molecule-wall collisions.
    fn handle_collision(&mut self) -> usize {
        let mut total_collisions = 0;
        // build a KD tree of all molecules
        let mut indexed_molecule = self.molecules.iter()
            .enumerate()
            .map(|(i, m)| IndexedMolecule { index: i, pos: m.pos, radius: m.mol_type.radius() })
            .collect::<Vec<_>>();
        let kd_tree = kd_tree::Node::new(&mut indexed_molecule, self.kd_tree_max_elem, kd_tree::KDNodeDivideBy::X);
        // handle collisions between molecules
        let num_molecules = self.molecules.len();
        for i in 0..num_molecules {
            let mol1 = self.molecules[i];
            let neighbors = kd_tree.query_circle(mol1.pos, mol1.mol_type.radius() * 2.0);
            for indexed_mol2 in neighbors {
                let j = indexed_mol2.index;
                let mol2 = &self.molecules[indexed_mol2.index];
                if i < j && Vec2::length(mol1.pos - mol2.pos) < mol1.mol_type.radius() + mol2.mol_type.radius() {
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
            let radius = m.mol_type.radius();
            if m.pos.x - radius < 0.0 {
                m.pos.x = 2.0 * radius - m.pos.x;
                m.vel.x = -m.vel.x;
                total_collisions += 1;
            }
            if m.pos.x + radius > self.width {
                m.pos.x = 2.0 * self.width - 2.0 * radius - m.pos.x;
                m.vel.x = -m.vel.x;
                total_collisions += 1;
            }
            if m.pos.y - radius < 0.0 {
                m.pos.y = 2.0 * radius - m.pos.y;
                m.vel.y = -m.vel.y;
                total_collisions += 1;
            }
            if m.pos.y + radius > self.height {
                m.pos.y = 2.0 * self.height - 2.0 * radius - m.pos.y;
                m.vel.y = -m.vel.y;
                total_collisions += 1;
            }
        }
        total_collisions
    }
}

impl<T: MoleculeType> Model for Model2<T> {
    type Type = T;
    type AdvanceReturnType = usize;

    fn construct(width: f32, height: f32, num_molecule: usize, constructor: impl FnMut(usize) ->  Molecule<Self::Type>) -> Self {
        let molecules = (0..num_molecule).map(constructor).collect();
        Self {
            molecules, width, height, kd_tree_max_elem: 10,
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
        for _ in 0..4 {
            // move all molecules along their speeds
            for m in self.molecules.iter_mut() {
                m.pos += m.vel * dt * 0.25;
            }
            // handle collisions
            self.handle_collision();
        }
        self.handle_collision()
    }
}