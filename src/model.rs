use eframe::egui::Vec2;

#[derive(Clone, Copy, Default, Debug)]
pub struct Molecule {
    pub pos: Vec2,
    pub vel: Vec2,
    pub radius: f32,
    pub mass: f32,
}

/// The abstract trait of models of gases in a perfect rectangular container.
pub trait Model {
    type AdvanceReturnType;
    /// The dimension (width, height) of the container
    fn dimension(&self) -> (f32, f32);
    /// Get the number of gas molecules in the container
    fn num_molecule(&self) -> usize;
    /// Advance the system by time dt (in seconds).
    fn advance(&mut self, dt: f32) -> Self::AdvanceReturnType;
    /// Get the summaries all molecules
    fn get_molecules(&self) -> impl Iterator<Item = Molecule>;
}