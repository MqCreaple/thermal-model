use eframe::egui::{Color32, Vec2, ViewportBuilder};
use model::MoleculeType;
use utils::color_interp;
use visualizer::PlotOptions;

mod utils;
mod model;
mod model_1;
mod model_2;
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
    fn mass(&self) -> f32 { 1.0 }

    fn radius(&self) -> f32 { 0.1 }

    fn color(&self, _pos: Vec2, vel: Vec2) -> eframe::egui::Color32 {
        let velocity = vel.length();
        let velocity_ratio = 1.0f32.min(velocity / MAX_VELOCITY);
        color_interp(COLOR_COLD, COLOR_HOT, velocity_ratio)
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::builder().filter_level(log::LevelFilter::Info).init();
    let model = model_2::Model2::<MoleculeType1>::new(120.0, 120.0, 15000, 10);
    // to change the effect of the visualizer, modify the options here
    let visualizer_options = visualizer::VisualizerOptions {
        plot_quantities: vec![
            (|m| m.vel.length(), "velocity magnitude"),
            (|m| m.vel.x, "velocity x"),
            (|m| m.vel.y, "velocity y"),
        ],
        plot_options: PlotOptions::Grid(50, 50),
    };

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
    let model = model_2::Model2::<MoleculeType1>::new(12.0, 12.0, 15000, 10);
    // to change the effect of the visualizer, modify the options here
    let visualizer_options = visualizer::VisualizerOptions {
        enable_plot_v_magnitude: true,
        enable_plot_v_x: true,
        enable_plot_v_y: true,
        plot_options: PlotOptions::Grid(50, 50),
        display_max_velocity: 1.414,
    };

    eframe::WebLogger::init(log::LevelFilter::Warn).unwrap();
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
