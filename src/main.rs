use eframe::egui::ViewportBuilder;
use visualizer::PlotOptions;

mod utils;
mod model;
mod model_1;
mod model_2;
mod visualizer;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model = model_2::Model2::new(80.0, 80.0, 0.1, 1.0, 8000, 10);
    let mut native_options = eframe::NativeOptions::default();
    native_options.viewport = ViewportBuilder::default().with_fullscreen(true);

    // to change the effect of the visualizer, modify the options here
    let visualizer_options = visualizer::VisualizerOptions {
        enable_plot_v_magnitude: true,
        enable_plot_v_x: true,
        enable_plot_v_y: true,
        plot_options: PlotOptions::Grid(50, 50),
        display_max_velocity: 1.414,
    };

    eframe::run_native(
        "gas molecule visualizer",
        native_options,
        Box::new(|cc| Box::new(
            visualizer::Visualizer::new(model, 0.25, 20, visualizer_options, &cc)
        ))
    )?;
    Ok(())
}
