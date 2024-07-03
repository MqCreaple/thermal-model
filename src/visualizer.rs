use std::fmt::Debug;
use std::time::{Duration, Instant};

use crate::model::Model;
use crate::utils::color_interp;
use eframe::{self, CreationContext};
use eframe::egui::{self, CentralPanel, Label, Mesh, Painter, Pos2, Rect, SidePanel, Ui, Vec2};
use egui_plot::{Bar, BarChart, Plot};

pub enum PlotOptions {
    All,                  // draw all molecules
    Grid(usize, usize),   // draw a grid of the specified size. Use different shades to represent the number of molecules in each grid
}

/// Customizable options for the visualizer
pub struct VisualizerOptions {
    pub enable_plot_v_magnitude: bool,
    pub plot_options: PlotOptions,  // the threshold of the number of molecules to switch between painting all molecules and painting a grid
    pub enable_plot_v_x: bool,
    pub enable_plot_v_y: bool,
    pub display_max_velocity: f32,
}

// constants of visualizer
const COLOR_HOT: egui::Color32 = egui::Color32::from_rgb(228, 47, 47);
const COLOR_COLD: egui::Color32 = egui::Color32::from_rgb(43, 110, 197);
// const COLOR_HOT: egui::Color32 = egui::Color32::from_rgb(0, 0, 0);
// const COLOR_COLD: egui::Color32 = egui::Color32::from_rgb(255, 255, 255);

pub struct Visualizer<M: Model> {
    model: M,
    speed_ratio: f32,   // the ratio of the speed of the model to the actual speed
    bar_cnt: usize,     // number of bars in the summary histogram
    last_frame: Instant,
    options: VisualizerOptions,
}

impl<M: Model> Visualizer<M> {
    pub fn new(model: M, speed_ratio: f32, bar_cnt: usize, options: VisualizerOptions, cc: &CreationContext<'_>) -> Self {
        Self { model, speed_ratio, last_frame: Instant::now(), bar_cnt, options }
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
        Bar::new((i as f64 + 0.5) * bin_width as f64, *cnt as f64)
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
    A: Debug {
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        // advance the model
        let now = Instant::now();
        let time_diff = (now - self.last_frame).as_secs_f32();
        let ret = self.model.advance(time_diff * self.speed_ratio);
        self.last_frame = now;
        // paint
        egui::CentralPanel::default().show(ctx, |ui| {
            // resize the window
            let (width, height) = self.model.dimension();
            let min_pixels_per_point = f32::min(ui.available_width() / width, ui.available_height() / height);
            // paint the molecules
            SidePanel::left("molecules").exact_width(width * min_pixels_per_point).show(ctx, |ui| {
                let painter = Painter::new(ui.ctx().clone(), ui.layer_id(), ui.available_rect_before_wrap());
                match self.options.plot_options {
                    PlotOptions::All => {
                        // if the number of molecules are sufficiently small, paint them all
                        let molecules = self.model.get_molecules();
                        molecules.for_each(|m| {
                            let pos = m.pos * min_pixels_per_point;
                            let pos = Pos2::new(pos.x, pos.y);
                            let radius = m.radius * min_pixels_per_point;
                            let velocity=  m.vel.length();
                            let velocity_ratio = 1.0f32.min(velocity / self.options.display_max_velocity);
                            painter.circle(pos, radius, color_interp(COLOR_COLD, COLOR_HOT, velocity_ratio), egui::Stroke::default());
                        });
                    },
                    PlotOptions::Grid(grid_x_count, grid_y_count) => {
                        // draw a grid and count the number of molecules in each grid
                        let mut grid = vec![0u32; grid_x_count * grid_y_count];
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
                                let multiplier = 1.0 - f32::min(1.0, (count as usize * grid_x_count * grid_y_count) as f32 / self.model.num_molecule() as f32);
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
            });
            CentralPanel::default().show(ctx, |ui| {
                ui.vertical(|ui| {
                    // show graphics info
                    ui.add(Label::new(format!("FPS: {}", (1.0 / time_diff) as i32)));
                    ui.add(Label::new(format!("Last return value from advance: {:?}", ret)));
                    // show summaries of the gas molecules
                    if self.options.enable_plot_v_magnitude {
                        // plot the histogram of the velocity magnitudes
                        let velocities = self.model.get_molecules().map(|m| m.vel.length()).collect::<Vec<_>>();
                        plot_bar(ui, &velocities, self.bar_cnt, "Velocity Histogram", "velocity (magnitude)", "count");
                    }
                    // plot the histogram of x components of the velocities
                    if self.options.enable_plot_v_x {
                        let x_velocities = self.model.get_molecules().map(|m| m.vel.x).collect::<Vec<_>>();
                        plot_bar(ui, &x_velocities, self.bar_cnt, "X Velocity Histogram", "velocity (x)", "count");
                    }
                    // plot the histogram of y components of the velocities
                    if self.options.enable_plot_v_y {
                        let y_velocities = self.model.get_molecules().map(|m| m.vel.y).collect::<Vec<_>>();
                        plot_bar(ui, &y_velocities, self.bar_cnt, "Y Velocity Histogram", "velocity (y)", "count");
                    }
                });
            });
        });
        // request repaint to ensure a minimum update frequency
        ctx.request_repaint_after(Duration::from_secs_f32(0.0f32.max(1.0 / 120.0 - time_diff)));
    }
}