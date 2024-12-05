use eframe::egui::Vec2;

pub trait MoleculeType: Clone + Copy {
    /// The maximum possible radius between two molecules.
    ///
    /// For example, if there are 3 types of molecules, each with radius 1.0, 2.0, and
    /// 3.0, then `MAX_RADIUS` should be 3.0.
    ///
    /// This is used in collision detection.
    const MAX_RADIUS: f32;
    fn mass(&self) -> f32;
    fn radius(&self) -> f32;

    /// The molecule's rotational inertia.
    fn inertia(&self) -> f32 {
        self.mass() * self.radius() * self.radius() / 2.0
    }
}

/// A simple molecule type.
///
/// Records the position, velocity, and angle of a molecule.
#[derive(Clone, Copy, Default, Debug)]
pub struct Molecule<T: MoleculeType> {
    pub pos: Vec2,
    pub vel: Vec2,

    /// The angle and angular velocity of the molecule.
    ///
    /// For single atom gases, this field is always `None`.
    pub orient: Option<(f32, f32)>,

    pub mol_type: T,
}

impl<T: MoleculeType> Molecule<T> {
    /// Advance the position and angle of the molecule by time `dt`.
    pub fn advance_pos_angle(&mut self, dt: f32) {
        use std::f32::consts::PI;

        if let Some((angle, angular_vel)) = &mut self.orient {
            self.pos += self.vel * dt;
            *angle += *angular_vel * dt;
            if *angle > 2.0 * PI {
                *angle -= 2.0 * PI;
            } else if *angle < 0.0 {
                *angle += 2.0 * PI;
            }
        } else {
            self.pos += self.vel * dt;
        }
    }
}

/// The abstract trait of models of gases in a perfect rectangular container.
pub trait Model {
    type Type: MoleculeType;
    type AdvanceReturnType;
    /// Create a new model with the given dimension and number of gas molecules.
    fn construct(
        width: f32,
        height: f32,
        num_molecule: usize,
        constructor: impl FnMut(usize) -> Molecule<Self::Type>,
    ) -> Self;
    /// The dimension (width, height) of the container
    fn dimension(&self) -> (f32, f32);
    /// Get the number of gas molecules in the container
    fn num_molecule(&self) -> usize;
    /// Advance the system by time dt (in seconds).
    fn advance(&mut self, dt: f32) -> Self::AdvanceReturnType;
    /// Get the summaries all molecules
    fn get_molecules(&self) -> impl Iterator<Item = Molecule<Self::Type>>;
}
