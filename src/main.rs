use eframe::egui::{Color32, Vec2, ViewportBuilder};
use model::{MoleculeType, Model};
use utils::color_interp;
use visualizer::PlotOptions;

mod utils;
mod model;
mod model_1;
mod model_2;
mod model_3;
mod visualizer;

#[derive(Clone, Copy, Default, Debug)]
struct MoleculeType1;

// molecule color constants
const COLOR_HOT: Color32 = Color32::from_rgb(228, 47, 47);
const COLOR_COLD: Color32 = Color32::from_rgb(43, 110, 197);
const MAX_VELOCITY: f32 = 1.414;
// const COLOR_HOT: Color32 = Color32::from_rgb(0, 0, 0);
// const COLOR_COLD: Color32 = Color32::from_rgb(255, 255, 255);

impl MoleculeType for MoleculeType1 {
    const MAX_RADIUS_BETWEEN_MOLECULES: f32 = 0.2;

    fn mass(&self) -> f32 { 1.0 }

    fn radius(&self) -> f32 { 0.1 }
}

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::builder().filter_level(log::LevelFilter::Info).init();
    let model = model_3::Model3::<MoleculeType1, 1000>::new(150.0, 150.0, 20000, 10);
    // to change the effect of the visualizer, modify the options here
    let visualizer_options = visualizer::VisualizerOptions::<model_3::Model3::<MoleculeType1, 1000>, f32> {
        plot_quantities: vec![
            (|m| m.vel.length(), "velocity magnitude"),
            (|m| m.vel.x, "velocity x"),
            (|m| m.vel.y, "velocity y"),
        ],
        state_quantities: vec![
            (|model| Some(model.total_energy()), "total energy"),
        ],
        plot_options: PlotOptions::Grid(50, 50),
        molecule_color: |m| {
            let velocity = m.vel.length();
            let velocity_ratio = 1.0f32.min(velocity / MAX_VELOCITY);
            color_interp(COLOR_COLD, COLOR_HOT, velocity_ratio)
        },
        grid_quantity: |_m| 1.0,   // Meaning that for every molecule in a grid, the grid will add 1 to its counter.
        grid_color: |count, model, x_count, y_count| {
            let multiplier = f32::min(1.0, 0.5 * (count as usize * x_count * y_count) as f32 / model.num_molecule() as f32);
            Color32::from_rgb((255.0 * multiplier) as u8, (255.0 * multiplier) as u8, (255.0 * multiplier) as u8)
        }
    };

    // start the app
    let mut native_options = eframe::NativeOptions::default();
    native_options.viewport = ViewportBuilder::default().with_fullscreen(true);
    eframe::run_native(
        "gas molecule visualizer",
        native_options,
        Box::new(|cc| Ok(Box::new(
            visualizer::Visualizer::new(model, 0.25, 20, visualizer_options, &cc)
        ))),
    )?;
    Ok(())
}

#[cfg(target_arch = "wasm32")]
fn main() {
    let model = model_2::Model2::<MoleculeType1>::new(120.0, 120.0, 15000, 10);
    // to change the effect of the visualizer, modify the options here
    let visualizer_options = visualizer::VisualizerOptions::<model_2::Model2<MoleculeType1>, f32> {
        plot_quantities: vec![
            (|m| m.vel.length(), "velocity magnitude"),
            (|m| m.vel.x, "velocity x"),
            (|m| m.vel.y, "velocity y"),
        ],
        state_quantities: vec![
            (|model| Some(model.total_energy()), "total energy"),
        ],
        plot_options: PlotOptions::Grid(50, 50),
        molecule_color: |m| {
            let velocity = m.vel.length();
            let velocity_ratio = 1.0f32.min(velocity / MAX_VELOCITY);
            color_interp(COLOR_COLD, COLOR_HOT, velocity_ratio)
        },
        grid_quantity: |_m| 1.0,   // Meaning that for every molecule in a grid, the grid will add 1 to its counter.
        grid_color: |count, model, x_count, y_count| {
            let multiplier = f32::min(1.0, 0.5 * (count as usize * x_count * y_count) as f32 / model.num_molecule() as f32);
            Color32::from_rgb((255.0 * multiplier) as u8, (255.0 * multiplier) as u8, (255.0 * multiplier) as u8)
        }
    };

    eframe::WebLogger::init(log::LevelFilter::Warn).expect("Failed to initialize web logger");
    let web_option = eframe::WebOptions::default();
    wasm_bindgen_futures::spawn_local(async move {
        eframe::WebRunner::new().start(
            "gas molecule visualizer",
            web_option,
            Box::new(|cc| Ok(Box::new(
                visualizer::Visualizer::new(model, 0.25, 20, visualizer_options, &cc)
            ))),
        )
        .await
        .expect("Failed to run web app")
    });
}
