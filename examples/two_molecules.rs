use eframe;
use eframe::egui::{Color32, Vec2, ViewportBuilder};
use rand::Rng;
use thermal_model::*;

/// Define the molecule type for O2 and CO2
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MoleculeTypes {
    O2, CO2,
}

const O2_COLOR: Color32 = Color32::from_rgb(228, 47, 47);
const CO2_COLOR: Color32 = Color32::from_rgb(32, 32, 32);

impl MoleculeType for MoleculeTypes {
    const MAX_RADIUS_BETWEEN_MOLECULES: f32 = 0.230 * 2.0;

    fn mass(&self) -> f32 {
        match self {
            MoleculeTypes::O2 => 32.0,
            MoleculeTypes::CO2 => 44.0,
        }
    }

    fn radius(&self) -> f32 {
        match self {
            MoleculeTypes::O2 => 0.152,
            MoleculeTypes::CO2 => 0.230,
        }
    }

    fn color(&self, _pos: Vec2, _vel: Vec2) -> eframe::egui::Color32 {
        match self {
            MoleculeTypes::O2 => O2_COLOR,
            MoleculeTypes::CO2 => CO2_COLOR,
        }
    
    }
}

const WIDTH: f32 = 150.0;
const HEIGHT: f32 = 150.0;
const NUM_MOLS: usize = 15000;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::builder().filter_level(log::LevelFilter::Info).init();
    let mut rng = rand::thread_rng();
    // generate some O2 and CO2 molecules
    let o2_mol_count = (NUM_MOLS as f32 * MoleculeTypes::O2.mass() / (MoleculeTypes::O2.mass() + MoleculeTypes::CO2.mass())) as usize;
    let model = Model2::<MoleculeTypes>::construct(WIDTH, HEIGHT, NUM_MOLS, |i| {
        if i < o2_mol_count {
            Molecule {
                pos: Vec2::new(rng.gen_range(0.0..WIDTH), rng.gen_range(0.0..HEIGHT)),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                mol_type: MoleculeTypes::O2,
            }
        } else {
            Molecule {
                pos: Vec2::new(rng.gen_range(0.0..WIDTH), rng.gen_range(0.0..HEIGHT)),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                mol_type: MoleculeTypes::CO2,
            }
        }
    });
    let mut native_options = eframe::NativeOptions::default();
    native_options.viewport = ViewportBuilder::default().with_fullscreen(true);
    eframe::run_native(
        "gas molecule visualizer",
        native_options,
        Box::new(|cc| Ok(Box::new(
            Visualizer::new(model, 0.25, 20, VisualizerOptions {
                plot_quantities: vec![
                    (|m| match m.mol_type {
                        MoleculeTypes::O2 => 0.5 * m.mol_type.mass() * m.vel.length_sq(),
                        MoleculeTypes::CO2 => f32::NAN,
                    }, "O2 kinetic energy"),
                    (|m| match m.mol_type {
                            MoleculeTypes::CO2 => 0.5 * m.mol_type.mass() * m.vel.length_sq(),
                            MoleculeTypes::O2 => f32::NAN,
                    }, "CO2 kinetic energy"),
                ],
                state_quantities: vec![
                    (|model| Some(model.total_energy() / model.num_molecule() as f32), "average kinetic energy"),
                    (|model| {
                        Some(
                            model.get_molecules()
                                .filter(|m| m.mol_type == MoleculeTypes::O2)
                                .fold(0.0, |acc, m| acc + 0.5 * m.mol_type.mass() * m.vel.length_sq())
                            / model.get_molecules().filter(|m| m.mol_type == MoleculeTypes::O2).count() as f32
                        )
                    }, "average energy of O2"),
                    (|model| {
                        Some(
                            model.get_molecules()
                                .filter(|m| m.mol_type == MoleculeTypes::CO2)
                                .fold(0.0, |acc, m| acc + 0.5 * m.mol_type.mass() * m.vel.length_sq())
                            / model.get_molecules().filter(|m| m.mol_type == MoleculeTypes::CO2).count() as f32
                        )
                    }, "average energy of CO2"),
                ],
                plot_options: PlotOptions::Grid(50, 50),
            }, &cc)
        ))),
    )?;
    Ok(())
}

#[cfg(target_arch = "wasm32")]
fn main() {
    let mut rng = rand::thread_rng();
    let o2_mol_count = (NUM_MOLS as f32 * MoleculeTypes::O2.mass() / (MoleculeTypes::O2.mass() + MoleculeTypes::CO2.mass())) as usize;
    let model = Model2::<MoleculeTypes>::construct(WIDTH, HEIGHT, NUM_MOLS, |i| {
        if i < o2_mol_count {
            Molecule {
                pos: Vec2::new(rng.gen_range(0.0..WIDTH), rng.gen_range(0.0..HEIGHT)),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                mol_type: MoleculeTypes::O2,
            }
        } else {
            Molecule {
                pos: Vec2::new(rng.gen_range(0.0..WIDTH), rng.gen_range(0.0..HEIGHT)),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                mol_type: MoleculeTypes::CO2,
            }
        }
    });
    eframe::WebLogger::init(log::LevelFilter::Warn).unwrap();
    let web_option = eframe::WebOptions::default();
    wasm_bindgen_futures::spawn_local(async move {
        eframe::WebRunner::new().start(
            "gas molecule visualizer",
            web_option,
            Box::new(|cc| Ok(Box::new(
                Visualizer::new(model, 0.25, 20, VisualizerOptions {
                    plot_quantities: vec![
                        (|m| match m.mol_type {
                            MoleculeTypes::O2 => 0.5 * m.mol_type.mass() * m.vel.length_sq(),
                            MoleculeTypes::CO2 => f32::NAN,
                        }, "O2 kinetic energy"),
                        (|m| match m.mol_type {
                                MoleculeTypes::CO2 => 0.5 * m.mol_type.mass() * m.vel.length_sq(),
                                MoleculeTypes::O2 => f32::NAN,
                        }, "CO2 kinetic energy"),
                    ],
                    state_quantities: vec![
                        (|model| Some(model.total_energy() / model.num_molecule() as f32), "average kinetic energy"),
                        (|model| {
                            Some(
                                model.get_molecules()
                                    .filter(|m| m.mol_type == MoleculeTypes::O2)
                                    .fold(0.0, |acc, m| acc + 0.5 * m.mol_type.mass() * m.vel.length_sq())
                                / model.get_molecules().filter(|m| m.mol_type == MoleculeTypes::O2).count() as f32
                            )
                        }, "average energy of O2"),
                        (|model| {
                            Some(
                                model.get_molecules()
                                    .filter(|m| m.mol_type == MoleculeTypes::CO2)
                                    .fold(0.0, |acc, m| acc + 0.5 * m.mol_type.mass() * m.vel.length_sq())
                                / model.get_molecules().filter(|m| m.mol_type == MoleculeTypes::CO2).count() as f32
                            )
                        }, "average energy of CO2"),
                    ],
                    plot_options: PlotOptions::Grid(50, 50),
                }, &cc)
            ))),
        )
        .await
        .expect("Failed to run web app")
    });
}