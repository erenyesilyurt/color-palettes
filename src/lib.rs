use color::OklabColor;
use color::OklabResult;
use color::RGBColor;
use plotters::prelude::BitMapBackend;
use plotters::prelude::DrawingBackend;
use std::cmp;
use std::error::Error;

pub mod color;

type Palette = Vec<RGBColor>;
const MAX_PALETTE_SIZE: u32 = 256;

pub enum GeneratorMode {
    Auto,
    // TODO: add more modes (monochromatic, complementary etc.)
}

pub struct GeneratorConfig {
    pub mode: GeneratorMode,
    pub similarity: f32,
    pub max_iterations: u32,
}

impl Default for GeneratorConfig {
    fn default() -> Self {
        Self {
            mode: GeneratorMode::Auto,
            similarity: 0.6,
            max_iterations: 50,
        }
    }
}

// Generate and add a color to the palette using farthest point sampling.
// step_size is the distance of the generated color to the last color in the palette
fn generate_with_fps(colors: &mut Palette, step_size: f32, max_iterations: u32) {
    let last_color = colors.last().unwrap().to_oklab();

    // Generate a set of candidate colors by taking random steps from the last color in the palette
    const POOL_SIZE: usize = 8; // generate up to 8 candidates
    let mut candidates = [OklabColor::new(); POOL_SIZE];
    {
        // try to generate unclipped candidates
        let mut num_candidates = 032;
        for i in 0..candidates.len() {
            let step_result = last_color.random_step_sgrb(step_size, max_iterations);
            if let OklabResult::Ok(c) = step_result {
                candidates[i] = c;
                num_candidates += 1;
            }
        }

        // generate and copy a single candidate if we failed above
        if num_candidates == 0 {
            candidates[0] = last_color
                .random_step_sgrb(step_size, max_iterations)
                .unwrap();
            for i in 1..candidates.len() {
                candidates[i] = candidates[0];
            }
        }
    }

    // Pick the candidate color that is the furthest away from existing colors
    let mut min_distances = [f32::INFINITY; POOL_SIZE];
    for (i, candidate) in candidates.iter().enumerate() {
        for color in colors.iter() {
            let color = color.to_oklab();
            let dist = candidate.squared_distance(&color);
            if dist < min_distances[i] {
                min_distances[i] = dist;
            }
        }
    }

    let mut max_idx: usize = 0;
    let mut max_dist: f32 = min_distances[0];
    for i in 1..min_distances.len() {
        if min_distances[i] > max_dist {
            max_dist = min_distances[i];
            max_idx = i;
        }
    }
    colors.push(RGBColor::from_oklab(candidates[max_idx]));
}

fn distance_from_similarity(similarity: f32) -> f32 {
    let similarity = similarity.clamp(0.0, 1.0);

    const MIN_DIST: f32 = 0.01;
    const MAX_DIST: f32 = 0.26;

    (MAX_DIST - similarity * (MAX_DIST - MIN_DIST)).max(0.0)
}

pub fn generate_palette(
    num_colors: u32,
    config: &GeneratorConfig,
) -> Result<Palette, Box<dyn Error>> {
    if num_colors == 0 || num_colors > MAX_PALETTE_SIZE {
        return Err(format!(
            "Number of colors must be number between  [1-{}].",
            MAX_PALETTE_SIZE
        )
        .into());
    }

    match config.mode {
        GeneratorMode::Auto => return Ok(generate_auto(num_colors, config)),
    }
}

fn generate_auto(num_colors: u32, config: &GeneratorConfig) -> Palette {
    let mut palette: Palette = Palette::with_capacity(num_colors as usize);
    let first_color = OklabColor::random_in_sgrb(0.3, 0.6, 100);
    palette.push(RGBColor::from_oklab(first_color.unwrap()));

    let color_distance = distance_from_similarity(config.similarity);
    for _ in 1..num_colors {
        generate_with_fps(&mut palette, color_distance, config.max_iterations);
    }

    palette
}

pub fn save_palette(palette: &Palette, path: &str) {
    if palette.len() == 0 {
        panic!("palette is length 0!");
    }

    let num_colors = palette.len() as u32;
    let squares_per_row = 256;
    let num_rows = (num_colors as f32 / squares_per_row as f32).ceil() as u32;

    let square_size: u32 = 32;
    let canvas_width = square_size * cmp::min(squares_per_row, num_colors);
    let canvas_height = square_size * num_rows;
    let mut backend = BitMapBackend::new(path, (canvas_width, canvas_height));

    let square_size = square_size as i32;
    for (i, &color) in palette.iter().enumerate() {
        let i = i as i32;
        let r: i32 = i / squares_per_row as i32;
        let c: i32 = i % squares_per_row as i32;

        let upper_left: (i32, i32) = (c * square_size, r * square_size);
        let bottom_right: (i32, i32) = (upper_left.0 + square_size, upper_left.1 + square_size);
        backend
            .draw_rect(
                upper_left,
                bottom_right,
                &plotters::style::RGBColor(color.r, color.g, color.b),
                true,
            )
            .unwrap();
    }
}
