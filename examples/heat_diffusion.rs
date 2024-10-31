use eframe::egui::ecolor::rgb_from_hsv;
use eframe::egui::{Color32, Vec2, ViewportBuilder};
use rand::Rng;
use std::f32;
use thermal_model::MoleculeType;
use thermal_model::{color_interp, linear_interp, sample_rand_velocity};
use thermal_model::{Model, Model3, Molecule};
use thermal_model::{PlotOptions, Visualizer, VisualizerOptions};

#[derive(Clone, Copy, Debug)]
struct MoleculeTypes;

impl MoleculeType for MoleculeTypes {
    const MAX_RADIUS: f32 = 0.1;

    fn mass(&self) -> f32 {
        1.0
    }

    fn radius(&self) -> f32 {
        0.1
    }
}

// model constants
const WIDTH: f32 = 300.0;
const HEIGHT: f32 = 300.0;
const NUM_MOLECULES: usize = 70000;
// visualizer constants
const COLOR_HOT: Color32 = Color32::from_rgb(228, 47, 47);
const COLOR_COLD: Color32 = Color32::from_rgb(43, 110, 197);
const MAX_VELOCITY: f32 = 1.414;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();
    let mut rng = rand::thread_rng();
    let model = Model3::<MoleculeTypes, 1>::construct(WIDTH, HEIGHT, NUM_MOLECULES, |i| {
        if i < NUM_MOLECULES / 2 {
            Molecule {
                pos: Vec2::new(
                    rng.gen_range(0.0..(WIDTH / 2.0)),
                    rng.gen_range(0.0..HEIGHT),
                ),
                vel: sample_rand_velocity(&mut rng, 1.0),
                mol_type: MoleculeTypes,
            }
        } else {
            Molecule {
                pos: Vec2::new(
                    rng.gen_range((WIDTH / 2.0)..WIDTH),
                    rng.gen_range(0.0..HEIGHT),
                ),
                vel: sample_rand_velocity(&mut rng, 0.4),
                mol_type: MoleculeTypes,
            }
        }
    });
    // to change the effect of the visualizer, modify the options here
    let visualizer_options = VisualizerOptions::<Model3<MoleculeTypes, 1>, Vec2> {
        plot_quantities: vec![(|m| m.vel.length(), "velocity magnitude")],
        state_quantities: vec![
            (|model| Some(model.total_energy()), "total energy"),
            (
                |model| {
                    let left = model.get_molecules().filter(|m| m.pos.x <= WIDTH / 2.0);
                    let mut count = 0;
                    let mut total_energy = 0.0;
                    for m in left {
                        count += 1;
                        total_energy += 0.5 * m.mol_type.mass() * m.vel.length_sq();
                    }
                    if count == 0 {
                        None
                    } else {
                        Some(total_energy / count as f32)
                    }
                },
                "average energy on the left half",
            ),
            (
                |model| {
                    let left = model.get_molecules().filter(|m| m.pos.x > WIDTH / 2.0);
                    let mut count = 0;
                    let mut total_energy = 0.0;
                    for m in left {
                        count += 1;
                        total_energy += 0.5 * m.mol_type.mass() * m.vel.length_sq();
                    }
                    if count == 0 {
                        None
                    } else {
                        Some(total_energy / count as f32)
                    }
                },
                "average energy on the right half",
            ),
        ],
        plot_options: PlotOptions::All,
        molecule_color: |m| {
            let velocity = m.vel.length();
            let velocity_ratio = 1.0f32.min(velocity / MAX_VELOCITY);
            color_interp(COLOR_COLD, COLOR_HOT, velocity_ratio)
        },
        grid_quantity: |m| Vec2::new(1.0, m.vel.length()),
        grid_color: |quantity, model, x_count, y_count| {
            let count = quantity.x;
            let total_velocity = quantity.y;
            let value = f32::min(
                1.0,
                0.5 * (count as usize * x_count * y_count) as f32 / model.num_molecule() as f32,
            );
            let avg_velocity = total_velocity / count;
            let velocity_ratio = 1.0f32.min(avg_velocity / MAX_VELOCITY);
            let hue = linear_interp(-1.0 / 3.0, 0.0, velocity_ratio);
            let saturation = 1.0 - 2.0 * velocity_ratio * (1.0 - velocity_ratio);
            let rgb = rgb_from_hsv((hue, saturation, value));
            Color32::from_rgb(
                (rgb[0] * 255.0) as u8,
                (rgb[1] * 255.0) as u8,
                (rgb[2] * 255.0) as u8,
            )
        },
    };

    // start the app
    let mut native_options = eframe::NativeOptions::default();
    native_options.viewport = ViewportBuilder::default().with_fullscreen(true);
    eframe::run_native(
        "gas molecule visualizer",
        native_options,
        Box::new(|cc| {
            Ok(Box::new(Visualizer::new(
                model,
                0.25,
                20,
                visualizer_options,
                &cc,
            )))
        }),
    )?;
    Ok(())
}
