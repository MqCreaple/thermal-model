use eframe::egui::{Color32, Vec2, ViewportBuilder};
use rand::Rng;
use thermal_model::MoleculeType;
use thermal_model::{color_interp, sample_rand_velocity};
use thermal_model::{Model, Model3, Molecule};
use thermal_model::{PlotOptions, Visualizer, VisualizerOptions};

#[derive(Clone, Copy, Debug)]
enum MoleculeTypes {
    Red,
    Green,
}

impl MoleculeType for MoleculeTypes {
    const MAX_RADIUS: f32 = 0.1;

    fn mass(&self) -> f32 {
        1.0
    }

    fn radius(&self) -> f32 {
        0.1
    }
}

const WIDTH: f32 = 200.0;
const HEIGHT: f32 = 200.0;
const NUM_MOLECULES: usize = 35000;
const RED_COLOR: Color32 = Color32::from_rgb(250, 10, 10);
const GREEN_COLOR: Color32 = Color32::from_rgb(0, 175, 20);

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
                vel: sample_rand_velocity(&mut rng, 0.64),
                orient: None,
                mol_type: MoleculeTypes::Red,
            }
        } else {
            Molecule {
                pos: Vec2::new(
                    rng.gen_range((WIDTH / 2.0)..WIDTH),
                    rng.gen_range(0.0..HEIGHT),
                ),
                vel: sample_rand_velocity(&mut rng, 0.64),
                orient: None,
                mol_type: MoleculeTypes::Green,
            }
        }
    });
    // to change the effect of the visualizer, modify the options here
    let visualizer_options = VisualizerOptions::<Model3<MoleculeTypes, 1>, Vec2> {
        plot_quantities: vec![
            (|m| m.vel.length(), "velocity magnitude"),
            (|m| m.vel.x, "velocity x"),
            (|m| m.vel.y, "velocity y"),
        ],
        state_quantities: vec![(|model| Some(model.translational_ke()), "total energy")],
        plot_options: PlotOptions::All,
        molecule_color: |m| match m.mol_type {
            MoleculeTypes::Red => RED_COLOR,
            MoleculeTypes::Green => GREEN_COLOR,
        },
        grid_quantity: |m| match m.mol_type {
            MoleculeTypes::Red => Vec2::new(1.0, 0.0),
            MoleculeTypes::Green => Vec2::new(0.0, 1.0),
        }, // Meaning that for every molecule in a grid, the grid will add 1 to its counter.
        grid_color: |quantity, _, _, _| {
            let red_cnt = quantity.x;
            let green_cnt = quantity.y;
            color_interp(RED_COLOR, GREEN_COLOR, red_cnt / (red_cnt + green_cnt))
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
