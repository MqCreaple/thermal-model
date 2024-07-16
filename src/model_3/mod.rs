use std::collections::BTreeSet;
use std::ops::{Add, Div};
use std::sync::{Arc, Mutex};
use std::{array, thread};

pub mod thread_sync;
use ordered_float::NotNan;
use thread_sync::ThreadSync;

use crate::model::{Model, Molecule, MoleculeType};
use crate::model_2::{kd_tree, IndexedMolecule};
use crate::utils;

use eframe::egui::Vec2;
use rand::Rng;
use ringbuffer::{ConstGenericRingBuffer, RingBuffer};

#[derive(PartialEq, Debug, Default, Clone, Copy)]
pub struct PressureEntry {
    pub x_neg: f32,
    pub x_pos: f32,
    pub y_neg: f32,
    pub y_pos: f32,
}

impl Add for PressureEntry {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x_neg: self.x_neg + rhs.x_neg,
            x_pos: self.x_pos + rhs.x_pos,
            y_neg: self.y_neg + rhs.y_neg,
            y_pos: self.y_pos + rhs.y_pos,
        }
    }

}

impl Div<f32> for PressureEntry {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        Self {
            x_neg: self.x_neg / rhs,
            x_pos: self.x_pos / rhs,
            y_neg: self.y_neg / rhs,
            y_pos: self.y_pos / rhs,
        }
    }

}

/// A wrapper type for `*mut Molecule<T>` that implements `Send`.
struct MoleculeMutPtr<T: MoleculeType>(*mut Molecule<T>);

unsafe impl<T: MoleculeType> Send for MoleculeMutPtr<T> {}

/// The type of entries stored in boundary sets
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SetKey(usize, NotNan<f32>);

impl PartialOrd for SetKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.1.partial_cmp(&other.1)
    }
}

impl Ord for SetKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.1.cmp(&other.1)
    }
}

/// The third model that records all the thermodynamic properties of the gas.
/// 
/// It also uses multithreading to improve the performance of the simulation.
/// 
/// `N` is the number of frames to record.
pub struct Model3<T: MoleculeType, const N: usize> {
    molecules: Vec<Molecule<T>>,
    width: f32,
    height: f32,
    unit_width: f32,    // width / THREAD_WIDTH
    unit_height: f32,   // height / THREAD_HEIGHT
    kd_tree_max_elem: usize,
    /// The ring buffer of pressures on the four walls
    pub(crate) pressures: ConstGenericRingBuffer<PressureEntry, N>,
}

pub const THREAD_WIDTH: usize = 4;
pub const THREAD_HEIGHT: usize = 4;

impl<T: MoleculeType + Send + Sync, const N: usize> Model3<T, N> {

    fn handle_collision(&mut self, dt: f32) -> PressureEntry {
        // divide the space to THREAD_WIDTH times THREAD_HEIGHT grids and
        // assign each grid a thread
        #[allow(const_evaluatable_unchecked)]
        thread::scope(|s| {
            let thread_sync = ThreadSync::new((THREAD_WIDTH * THREAD_HEIGHT) as u32);
            let _panic_hook_guard = thread_sync.set_panic_hook();
            let num_molecules = self.molecules.len();
            let kd_tree_max_elem = self.kd_tree_max_elem;

            let horizontal_bound_mols: [[Arc<Mutex<BTreeSet<SetKey>>>; THREAD_HEIGHT]; THREAD_WIDTH + 1] =
                array::from_fn(|_| array::from_fn(|_| Arc::new(Mutex::new(BTreeSet::new()))));
            let vertical_bound_mols: [[Arc<Mutex<BTreeSet<SetKey>>>; THREAD_HEIGHT + 1]; THREAD_WIDTH] =
                array::from_fn(|_| array::from_fn(|_| Arc::new(Mutex::new(BTreeSet::new()))));

            for thread_idx in 0..THREAD_WIDTH * THREAD_HEIGHT {
                let molecules = MoleculeMutPtr(self.molecules.as_mut_ptr());
                let thread_sync_clone = thread_sync.clone();
                let grid_i = thread_idx % THREAD_WIDTH;
                let grid_j = thread_idx / THREAD_WIDTH;
                let grid_x_range = (grid_i as f32 * self.unit_width)..((grid_i + 1) as f32 * self.unit_width);
                let grid_y_range = (grid_j as f32 * self.unit_height)..((grid_j + 1) as f32 * self.unit_height);
                let left_bound_mols = horizontal_bound_mols[grid_i][grid_j].clone();
                let right_bound_mols = horizontal_bound_mols[grid_i + 1][grid_j].clone();
                let up_bound_mols = vertical_bound_mols[grid_i][grid_j].clone();
                let down_bound_mols = vertical_bound_mols[grid_i][grid_j + 1].clone();
                s.spawn(move || {
                    let _ = &molecules;
                    let mut indexed_molecules = (0..num_molecules)
                        .map(|i| (i, unsafe { molecules.0.add(i) }))
                        .filter(|&(_i, m)| unsafe {
                                grid_x_range.contains(&(*m).pos.x)
                                && grid_y_range.contains(&(*m).pos.y)
                        })
                        .map(|(i, m)| unsafe {
                            IndexedMolecule { index: i, pos: (*m).pos }
                        })
                        .collect::<Vec<_>>();
                    let index_set = indexed_molecules.iter().map(|m| m.index).collect::<Vec<_>>();
                    // TODO URGENT: sometimes deadlocks happen in the code. Need to debug.
                    thread_sync_clone.sync().unwrap();
                    log::debug!("First sync:\t{}", thread_idx);
                    // This operation may not look safe, but since before the sync statement,
                    // every thread saves a different set of indices to check collision.
                    // In other words, let M be the set of molecules, then thread i is selecting
                    // a subset of indices M_i, such that union of all M_i is M and intersection
                    // of any M_i with M_j is null.
                    // Therefore, no data race will happen, even if the M_i is different for each
                    // frame.
                    let kd_tree = kd_tree::Node::new(&mut indexed_molecules, kd_tree_max_elem, kd_tree::KDNodeDivideBy::X);
                    // process collisions between molecules within the grid
                    for i in &index_set {
                        let i = *i;
                        let mol1 = unsafe { molecules.0.add(i).as_mut() };
                        if let Some(mol1) = mol1 {
                            let radius1 = mol1.mol_type.radius();
                            let neighbors = kd_tree.query_circle(mol1.pos, T::MAX_RADIUS_BETWEEN_MOLECULES).unwrap_or(Vec::new());
                            for mol2 in neighbors {
                                let j = mol2.index;
                                let mol2 = unsafe { molecules.0.add(j).as_mut() };
                                if let Some(mol2) = mol2 {
                                    let radius2 = mol2.mol_type.radius();
                                    if i < j && Vec2::length(mol1.pos - mol2.pos) < radius1 + radius2 {
                                        utils::collision_2_molecules(mol1, mol2);
                                    }
                                } else {
                                    continue;
                                }
                            }
                        } else {
                            continue;
                        }
                    }
                    // process collisions with the grid boundary
                    // add the molecules that are on the boundaries to `horizontal_bound_mols` and `vertical_bound_mols`
                    for i in &index_set {
                        let i = *i;
                        let m = unsafe { *molecules.0.add(i) };
                        let radius = m.mol_type.radius();
                        // TODO: mutex lock is time inefficient
                        if m.pos.x - radius < grid_x_range.start {
                            left_bound_mols.lock().unwrap().insert(SetKey(i, NotNan::new(m.pos.y).unwrap()));
                        }
                        if m.pos.x + radius > grid_x_range.end {
                            right_bound_mols.lock().unwrap().insert(SetKey(i, NotNan::new(m.pos.y).unwrap()));
                        }
                        if m.pos.y - radius < grid_y_range.start {
                            up_bound_mols.lock().unwrap().insert(SetKey(i, NotNan::new(m.pos.x).unwrap()));
                        }
                        if m.pos.y + radius > grid_y_range.end {
                            down_bound_mols.lock().unwrap().insert(SetKey(i, NotNan::new(m.pos.x).unwrap()));
                        }
                    }
                    thread_sync_clone.sync().unwrap();
                    log::debug!("Second sync:\t{}", thread_idx);

                    // each thread picks a boundary and process collision
                    if grid_i >= 1 {
                        // evaluate collisions on horizontal bound
                        let bound_mols = left_bound_mols.lock().unwrap();
                        for key1 in bound_mols.iter() {
                            let i = key1.0;
                            let mol1 = unsafe { molecules.0.add(i).as_mut().unwrap() };
                            let radius1 = mol1.mol_type.radius();
                            let up_bound = NotNan::new(mol1.pos.y - T::MAX_RADIUS_BETWEEN_MOLECULES).unwrap();
                            let down_bound = NotNan::new(mol1.pos.y + T::MAX_RADIUS_BETWEEN_MOLECULES).unwrap();
                            for key2 in bound_mols.range(SetKey(i, up_bound)..SetKey(i, down_bound)) {
                                let j = key2.0;
                                let mol2 = unsafe { molecules.0.add(j).as_mut().unwrap() };
                                let radius2 = mol2.mol_type.radius();
                                if Vec2::length(mol1.pos - mol2.pos) <= radius1 + radius2 {
                                    utils::collision_2_molecules(mol1, mol2);
                                }
                            }
                        }
                    }
                    if grid_j >= 1 {
                        let bound_mols = up_bound_mols.lock().unwrap();
                        for key1 in bound_mols.iter() {
                            let i = key1.0;
                            let mol1 = unsafe { molecules.0.add(i).as_mut().unwrap() };
                            let radius1 = mol1.mol_type.radius();
                            let left_bound = NotNan::new(mol1.pos.x - T::MAX_RADIUS_BETWEEN_MOLECULES).unwrap();
                            let right_bound = NotNan::new(mol1.pos.x + T::MAX_RADIUS_BETWEEN_MOLECULES).unwrap();
                            for key2 in bound_mols.range(SetKey(i, left_bound)..SetKey(i, right_bound)) {
                                let j = key2.0;
                                let mol2 = unsafe { molecules.0.add(j).as_mut().unwrap() };
                                let radius2 = mol2.mol_type.radius();
                                if Vec2::length(mol1.pos - mol2.pos) <= radius1 + radius2 {
                                    utils::collision_2_molecules(mol1, mol2);
                                }
                            }
                        }
                    }
                });
            }
        });

        // handle collisions with walls
        let mut pressure = PressureEntry::default();
        for m in self.molecules.iter_mut() {
            let radius = m.mol_type.radius();
            if m.pos.x - radius < 0.0 {
                m.pos.x = 2.0 * radius - m.pos.x;
                m.vel.x = -m.vel.x;
                pressure.x_neg += 2.0 * m.vel.x * m.mol_type.mass();
            }
            if m.pos.x + radius > self.width {
                m.pos.x = 2.0 * self.width - 2.0 * radius - m.pos.x;
                m.vel.x = -m.vel.x;
                pressure.x_pos -= 2.0 * m.vel.x * m.mol_type.mass();
            }
            if m.pos.y - radius < 0.0 {
                m.pos.y = 2.0 * radius - m.pos.y;
                m.vel.y = -m.vel.y;
                pressure.y_neg += 2.0 * m.vel.y * m.mol_type.mass();
            }
            if m.pos.y + radius > self.height {
                m.pos.y = 2.0 * self.height - 2.0 * radius - m.pos.y;
                m.vel.y = -m.vel.y;
                pressure.y_pos -= 2.0 * m.vel.y * m.mol_type.mass();
            }
        }
        pressure.x_neg /= self.width * dt;
        pressure.x_pos /= self.width * dt;
        pressure.y_neg /= self.height * dt;
        pressure.y_pos /= self.height * dt;
        pressure
    }

    pub fn volume(&self) -> f32 {
        self.width * self.height
    }

    pub fn total_energy(&self) -> f32 {
        self.molecules.iter().fold(0.0, |acc, m| {
            acc + 0.5 * m.mol_type.mass() * m.vel.length_sq()
        })
    }

    /// Get the pressure of the gas
    /// 
    /// If the pressure on the four walls differs by more than 5%, then return None.
    pub fn average_pressure(&self) -> Option<f32> {
        let pressure_entry = self.pressures.iter().fold(PressureEntry::default(), |acc, p| acc + *p) / (self.pressures.len() as f32);
        let pressure = (pressure_entry.x_neg + pressure_entry.x_pos + pressure_entry.y_neg + pressure_entry.y_pos) / 4.0;
        if (pressure_entry.x_neg - pressure).abs() / pressure > 0.05
            || (pressure_entry.x_pos - pressure).abs() / pressure > 0.05
            || (pressure_entry.y_neg - pressure).abs() / pressure > 0.05
            || (pressure_entry.y_pos - pressure).abs() / pressure > 0.05 {
            None
        } else {
            Some(pressure)
        }
    }

    pub fn pressure(&self) -> PressureEntry {
        self.pressures.iter().fold(PressureEntry::default(), |acc, p| acc + *p) / (N as f32)
    }
}

impl<T: MoleculeType + Default + Send + Sync, const N: usize> Model3<T, N> {
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
            unit_width: width / THREAD_WIDTH as f32,
            unit_height: height / THREAD_HEIGHT as f32,
            pressures: ConstGenericRingBuffer::new(),
        }
    }
}

impl<T: MoleculeType + Send + Sync, const N: usize> Model for Model3<T, N> {
    type Type = T;
    type AdvanceReturnType = usize;

    fn construct(width: f32, height: f32, num_molecule: usize, constructor: impl FnMut(usize) ->  Molecule<Self::Type>) -> Self {
        let molecules = (0..num_molecule).map(constructor).collect();
        Self {
            molecules, width, height, kd_tree_max_elem: 10,
            unit_width: width / THREAD_WIDTH as f32,
            unit_height: height / THREAD_HEIGHT as f32,
            pressures: ConstGenericRingBuffer::new(),
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
            let pressure = self.handle_collision(dt * 0.25);
            self.pressures.push(pressure);
        }
        0
    }
}