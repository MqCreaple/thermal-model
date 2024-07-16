use eframe::egui::{Color32, Vec2, ViewportBuilder};
use thermal_model::{color_interp, MoleculeType};
use thermal_model::{Model3, Visualizer, VisualizerOptions};

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

    fn color(&self, _pos: Vec2, vel: Vec2) -> eframe::egui::Color32 {
        let velocity = vel.length();
        let velocity_ratio = 1.0f32.min(velocity / MAX_VELOCITY);
        color_interp(COLOR_COLD, COLOR_HOT, velocity_ratio)
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::builder().filter_level(log::LevelFilter::Info).init();
    let model = Model3::<MoleculeType1, 10000>::new(300.0, 300.0, 80000, 10);
    // to change the effect of the visualizer, modify the options here
    let visualizer_options = VisualizerOptions::<Model3<MoleculeType1, 10000>> {
        plot_quantities: vec![
            (|m| m.vel.length(), "velocity magnitude"),
            (|m| m.vel.x, "velocity x"),
            (|m| m.vel.y, "velocity y"),
        ],
        state_quantities: vec![
            (|model| Some(model.total_energy()), "total energy"),
            (|model| model.average_pressure(), "average pressure"),
            (|model| Some(model.pressure().x_pos), "pressure on positive X side"),
            (|model| Some(model.pressure().x_neg), "pressure on negative X side"),
            (|model| Some(model.pressure().y_pos), "pressure on positive Y side"),
            (|model| Some(model.pressure().y_neg), "pressure on negative Y side"),
        ],
        plot_options: thermal_model::PlotOptions::Grid(50, 50),
    };

    let mut native_options = eframe::NativeOptions::default();
    native_options.viewport = ViewportBuilder::default().with_fullscreen(true);
    eframe::run_native(
        "gas molecule visualizer",
        native_options,
        Box::new(|cc| Ok(Box::new(
            Visualizer::new(model, 0.25, 20, visualizer_options, &cc)
        ))),
    )?;
    Ok(())
}

#[cfg(target_arch = "wasm32")]
fn main() {
    let model = Model3::<MoleculeType1>::new(120.0, 120.0, 15000, 10);
    // to change the effect of the visualizer, modify the options here
    let visualizer_options = VisualizerOptions::<Model3<MoleculeType1, 10000>> {
        plot_quantities: vec![
            (|m| m.vel.length(), "velocity magnitude"),
            (|m| m.vel.x, "velocity x"),
            (|m| m.vel.y, "velocity y"),
        ],
        state_quantities: vec![
            (|model| Some(model.total_energy()), "total energy"),
            (|model| model.average_pressure(), "average pressure"),
            (|model| Some(model.raw_pressure().x_pos), "pressure on positive X side"),
            (|model| Some(model.raw_pressure().x_neg), "pressure on negative X side"),
            (|model| Some(model.raw_pressure().y_pos), "pressure on positive Y side"),
            (|model| Some(model.raw_pressure().y_neg), "pressure on negative Y side"),
        ],
        plot_options: thermal_model::PlotOptions::Grid(50, 50),
    };

    eframe::WebLogger::init(log::LevelFilter::Warn).unwrap();
    let web_option = eframe::WebOptions::default();
    wasm_bindgen_futures::spawn_local(async move {
        eframe::WebRunner::new().start(
            "gas molecule visualizer",
            web_option,
            Box::new(|cc| Ok(Box::new(
                Visualizer::new(model, 0.25, 20, visualizer_options, &cc)
            ))),
        )
        .await
        .expect("Failed to run web app")
    });
}
