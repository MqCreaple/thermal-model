use std::fmt::Debug;
use std::time::{Duration, Instant};

use crate::model::{Model, MoleculeType};
use crate::model::Molecule;
use chrono::{DateTime, Utc};
use eframe::{self, CreationContext};
use eframe::egui::{self, CentralPanel, Color32, Grid, Label, Mesh, Painter, Pos2, Rect, SidePanel, Stroke, TopBottomPanel, Ui, Vec2};
use egui_plot::{Bar, BarChart, Plot};
use log::debug;
use ordered_float::NotNan;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlotOptions {
    All,                  // draw all molecules
    Grid(usize, usize),   // draw a grid of the specified size. Use different shades to represent the number of molecules in each grid
}

/// Customizable options for the visualizer
pub struct VisualizerOptions<M: Model> {
    /// The way to show the molecules
    pub plot_options: PlotOptions,

    /// The quantities to take a statistical summary of, along with each quantity's name
    pub plot_quantities: Vec<(fn(Molecule<M::Type>) -> f32, &'static str)>,

    /// The state quantities (e.g. pressure, temperature, entropy) to display.
    ///
    /// `None` represents a non-well-defined quantity
    pub state_quantities: Vec<(fn(&M) -> Option<f32>, &'static str)>,                
}

pub struct Visualizer<M: Model> {
    model: M,
    speed_ratio: f32,   // the ratio of the speed of the model to the actual speed
    bar_cnt: usize,     // number of bars in the summary histogram
    last_frame: DateTime<Utc>,
    options: VisualizerOptions<M>,
}

impl<M: Model> Visualizer<M> {
    pub const FPS: f32 = 120.0;

    pub fn new(model: M, speed_ratio: f32, bar_cnt: usize, options: VisualizerOptions<M>, cc: &CreationContext<'_>) -> Self {
        Self { model, speed_ratio, last_frame: Utc::now(), bar_cnt, options }
    }
}

fn plot_bar<'a>(ui: &mut Ui, data: &[f32], bar_cnt: usize, id: &'static str, x_label: &'static str, y_label: &'static str) {
    let min_v = data.iter().fold(f32::INFINITY, |acc, &v| acc.min(v));
    let max_v = data.iter().fold(0.0f32, |acc, &v| acc.max(v));
    let mut bins = vec![0u32; bar_cnt];
    let bin_width = (max_v - min_v) / bar_cnt as f32;
    data.iter().for_each(|&v| {
        let bin = f32::min((v - min_v) / bin_width, (bar_cnt - 1) as f32) as usize;
        bins[bin] += 1;
    });
    let bars = bins.iter().enumerate().map(|(i, cnt)|
        Bar::new(min_v as f64 + (i as f64 + 0.5) * bin_width as f64, *cnt as f64)
    ).collect();
    let bar_chart = BarChart::new(bars).width(bin_width as f64);
    Plot::new(id)
        .allow_drag(false).allow_scroll(false).allow_boxed_zoom(false).allow_zoom(false).allow_double_click_reset(false)
        .x_axis_label(x_label)
        .y_axis_label(y_label)
        .view_aspect(1.6)
        .show(ui, |plot_ui| plot_ui.bar_chart(bar_chart));
}

impl<M, A> eframe::App for Visualizer<M>
where
    M: Model<AdvanceReturnType = A>,
    A: Debug
{
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        // advance the model
        let now = Utc::now();
        let time_diff = now - self.last_frame;
        let time_diff = time_diff.num_seconds() as f32 + time_diff.subsec_nanos() as f32 / 1_000_000_000.0;
        let ret = self.model.advance(time_diff * self.speed_ratio);
        self.last_frame = now;
        // paint
        egui::CentralPanel::default().show(ctx, |ui| {
            // resize the window
            let (width, height) = self.model.dimension();
            let min_pixels_per_point = f32::min(ui.available_width() / width, ui.available_height() / height);
            // if the width is larger than the height, use a left-right layout
            // otherwise, use a top-bottom layout
            let use_horizontal_layout = ui.available_width() / width > ui.available_height() / height;
            // paint the molecules
            let paint_molecules = |ui: &mut Ui| {
                let painter = Painter::new(
                    ui.ctx().clone(),
                    ui.layer_id(),
                    Rect::from_min_size(Pos2::ZERO, Vec2::new(width * min_pixels_per_point, height * min_pixels_per_point))
                );
                match self.options.plot_options {
                    PlotOptions::All => {
                        // if the number of molecules are sufficiently small, paint them all
                        let molecules = self.model.get_molecules();
                        molecules.for_each(|m| {
                            let pos = m.pos * min_pixels_per_point;
                            let pos = Pos2::new(pos.x, pos.y);
                            let radius = m.mol_type.radius() * min_pixels_per_point;
                            painter.circle(pos, radius, m.mol_type.color(m.pos, m.vel), egui::Stroke::default());
                        });
                        // paint the boundary
                        painter.rect(
                            Rect::from_min_size(Pos2::ZERO, Vec2::new(width * min_pixels_per_point, height * min_pixels_per_point)),
                            0.0,
                            Color32::TRANSPARENT,
                            Stroke::new(1.0, Color32::BLACK)
                        );
                    },
                    PlotOptions::Grid(grid_x_count, grid_y_count) => {
                        // draw a grid and count the number of molecules in each grid
                        let grid_size = grid_x_count * grid_y_count;
                        let mut grid = vec![0u32; grid_size];
                        let molecules = self.model.get_molecules();
                        let grid_width = width / grid_x_count as f32;
                        let grid_height = height / grid_y_count as f32;
                        molecules.for_each(|m| {
                            let x = (m.pos.x / grid_width).floor() as usize;
                            let y = (m.pos.y / grid_height).floor() as usize;
                            if x < grid_x_count && y < grid_y_count {
                                grid[x * grid_y_count + y] += 1;
                            }
                        });
                        let mut mesh = Mesh::default();
                        let display_width = width * min_pixels_per_point;
                        let display_height = height * min_pixels_per_point;
                        for x in 0..grid_x_count {
                            for y in 0..grid_y_count {
                                let count = grid[x * grid_y_count + y];
                                let multiplier = f32::min(1.0, 0.5 * (count as usize * grid_size) as f32 / self.model.num_molecule() as f32);
                                let color = egui::Color32::from_rgb((255.0 * multiplier) as u8, (255.0 * multiplier) as u8, (255.0 * multiplier) as u8);
                                let x = x as f32 * display_width / grid_x_count as f32;
                                let y = y as f32 * display_height / grid_y_count as f32;
                                let width = display_width / grid_x_count as f32;
                                let height = display_height / grid_y_count as f32;
                                mesh.add_colored_rect(
                                    Rect::from_min_size(Pos2::new(x, y), Vec2::new(width, height)),
                                    color
                                );
                            }
                        }
                        painter.add(egui::Shape::mesh(mesh));
                    }
                }
            };

            // handle dynamic layout
            if use_horizontal_layout {
                SidePanel::left("left panel").exact_width(min_pixels_per_point * width).show(ctx, paint_molecules);
            } else {
                TopBottomPanel::top("top panel").exact_height(min_pixels_per_point * height).show(ctx, paint_molecules);
            }

            CentralPanel::default().show(ctx, |ui| {
                ui.vertical(|ui| {
                    // show graphics info
                    ui.add(Label::new(format!("FPS: {}", (1.0 / time_diff) as i32)));
                    ui.add(Label::new(format!("Number of molecules: {}", self.model.num_molecule())));
                    ui.add(Label::new(format!("Last number of unchecked collisions: {:?}", ret)));
                    
                    // show control widgets
                    let mut plot_all_molecules = self.options.plot_options == PlotOptions::All;
                    let (mut grid_x_count, mut grid_y_count) = match self.options.plot_options {
                        PlotOptions::Grid(x, y) => (x, y),
                        _ => (50, 50),
                    };
                    ui.group(|ui| {
                        ui.checkbox(&mut plot_all_molecules, "Plot all molecules");
                        ui.horizontal(|ui| {
                            ui.add_enabled(!plot_all_molecules, egui::Slider::new(&mut grid_x_count, 1..=100).text("Grid X Count"));
                            ui.add_enabled(!plot_all_molecules, egui::Slider::new(&mut grid_y_count, 1..=100).text("Grid Y Count"));
                        });
                        if !plot_all_molecules {
                            self.options.plot_options = PlotOptions::Grid(grid_x_count, grid_y_count);
                        } else {
                            self.options.plot_options = PlotOptions::All;
                        }
                    });
        
                    // show all state variables
                    for (state_fn, state_name) in &self.options.state_quantities {
                        let quantity = state_fn(&self.model);
                        let label = if let Some(quantity) = quantity {
                            format!("{}: {}", state_name, quantity)
                        } else {
                            format!("{}: Not well defined.", state_name)
                        };
                        ui.add(Label::new(label));
                    }
                    
                    // show all plots of quantites to plot
                    let mut graph_grid = Grid::new("graphs");
                    if use_horizontal_layout {
                        graph_grid = graph_grid.min_col_width(ui.available_width());
                    } else {
                        graph_grid = graph_grid.min_row_height(ui.available_height());
                    }
                    graph_grid.show(ui, |ui| {
                        for quantity in self.options.plot_quantities.iter() {
                            let (f, id) = quantity;
                            // NaN data represents absense. Filter out NaN data before plotting.
                            let data = self.model.get_molecules().map(f).filter(|data| !data.is_nan()).collect::<Vec<f32>>();
                            plot_bar(ui, &data, self.bar_cnt, id, id, "count");
                            if use_horizontal_layout { ui.end_row(); }
                        }
                    });
                });
            });
        });
        // request repaint to ensure a minimum update frequency
        ctx.request_repaint_after(Duration::from_secs_f32(0.0f32.max(1.0 / Self::FPS - time_diff)));
    }
}