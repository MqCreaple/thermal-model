use eframe;
use eframe::egui::{Color32, Vec2, ViewportBuilder};
use rand::Rng;
use thermal_model::*;

/// Define the molecule type for O2 and CO2
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MoleculeTypes {
    O2,
    CO2,
}

const O2_COLOR: Color32 = Color32::from_rgb(228, 47, 47);
const CO2_COLOR: Color32 = Color32::from_rgb(32, 32, 32);

impl MoleculeType for MoleculeTypes {
    const MAX_RADIUS: f32 = 0.230;

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
}

const WIDTH: f32 = 400.0;
const HEIGHT: f32 = 400.0;
const NUM_MOLS: usize = 80000;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use std::f32::consts::PI;

    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();
    let mut rng = rand::thread_rng();
    // generate some O2 and CO2 molecules
    let o2_mol_count = (NUM_MOLS as f32 * MoleculeTypes::O2.mass()
        / (MoleculeTypes::O2.mass() + MoleculeTypes::CO2.mass())) as usize;
    let model = Model3::<MoleculeTypes, 10000>::construct(WIDTH, HEIGHT, NUM_MOLS, |i| {
        if i < o2_mol_count {
            Molecule {
                pos: Vec2::new(rng.gen_range(0.0..WIDTH), rng.gen_range(0.0..HEIGHT)),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                mol_type: MoleculeTypes::O2,
                orient: Some((rng.gen_range(0.0..2.0 * PI), rng.gen_range(-1.0..1.0))),
            }
        } else {
            Molecule {
                pos: Vec2::new(rng.gen_range(0.0..WIDTH), rng.gen_range(0.0..HEIGHT)),
                vel: Vec2::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)),
                mol_type: MoleculeTypes::CO2,
                orient: Some((rng.gen_range(0.0..2.0 * PI), rng.gen_range(-1.0..1.0))),
            }
        }
    });

    let visualizer_options = VisualizerOptions::<Model3<MoleculeTypes, 10000>, f32> {
        plot_quantities: vec![
            (
                |m| match m.mol_type {
                    MoleculeTypes::O2 => 0.5 * m.mol_type.mass() * m.vel.length_sq(),
                    MoleculeTypes::CO2 => f32::NAN,
                },
                "O2 translational kinetic energy",
            ),
            (
                |m| match m.mol_type {
                    MoleculeTypes::CO2 => 0.5 * m.mol_type.mass() * m.vel.length_sq(),
                    MoleculeTypes::O2 => f32::NAN,
                },
                "CO2 translational kinetic energy",
            ),
        ],
        state_quantities: vec![
            (
                |model| Some(model.translational_ke() / model.num_molecule() as f32),
                "average translational kinetic energy",
            ),
            (
                |model| {
                    Some(
                        model
                            .get_molecules()
                            .filter(|m| m.mol_type == MoleculeTypes::O2)
                            .fold(0.0, |acc, m| {
                                acc + 0.5 * m.mol_type.mass() * m.vel.length_sq()
                            })
                            / model
                                .get_molecules()
                                .filter(|m| m.mol_type == MoleculeTypes::O2)
                                .count() as f32,
                    )
                },
                "average translational kinetic energy of O2",
            ),
            (
                |model| {
                    Some(
                        model
                            .get_molecules()
                            .filter(|m| m.mol_type == MoleculeTypes::CO2)
                            .fold(0.0, |acc, m| {
                                acc + 0.5 * m.mol_type.mass() * m.vel.length_sq()
                            })
                            / model
                                .get_molecules()
                                .filter(|m| m.mol_type == MoleculeTypes::CO2)
                                .count() as f32,
                    )
                },
                "average translational kinetic energy of CO2",
            ),
        ],
        plot_options: PlotOptions::All,
        molecule_color: |m| match m.mol_type {
            MoleculeTypes::CO2 => CO2_COLOR,
            MoleculeTypes::O2 => O2_COLOR,
        },
        grid_quantity: |_m| 1.0, // Meaning that for every molecule in a grid, the grid will add 1 to its counter.
        grid_color: |count, model, x_count, y_count| {
            let multiplier = f32::min(
                1.0,
                0.5 * (count as usize * x_count * y_count) as f32 / model.num_molecule() as f32,
            );
            Color32::from_rgb(
                (255.0 * multiplier) as u8,
                (255.0 * multiplier) as u8,
                (255.0 * multiplier) as u8,
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
