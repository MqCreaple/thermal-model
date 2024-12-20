mod model;
mod utils;

// Simple O(N^2) collision detection
mod model_1;
// KD tree for O(N log N) collision detection
mod model_2;
// KD tree + Multithreading
mod model_3;
mod visualizer;

pub use model::*;
pub use model_1::*;
pub use model_2::*;
pub use model_3::*;
pub use utils::*;
pub use visualizer::*;
