use eframe::egui::{Vec2, Color32};

use crate::model::Molecule;

pub fn collision_2_molecules(mol1: &mut Molecule, mol2: &mut Molecule) {
    // calculate velocities after collision
    let normal = (mol2.pos - mol1.pos).normalized();
    let m1 = mol1.mass;
    let m2 = mol2.mass;
    let v1 = mol1.vel.dot(normal);
    let v2 = mol2.vel.dot(normal);
    let v1_final = 2.0 * m2 * (v2 - v1) / (m1 + m2) * normal + mol1.vel;
    let v2_final = 2.0 * m1 * (v1 - v2) / (m1 + m2) * normal + mol2.vel;
    // update the two molecules' positions and velocities
    let update_dist = (mol1.radius + mol2.radius - Vec2::length(mol1.pos - mol2.pos)) / 2.0;
    mol1.vel = v1_final;
    mol2.vel = v2_final;
    mol1.pos -= normal * update_dist;
    mol2.pos += normal * update_dist;
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