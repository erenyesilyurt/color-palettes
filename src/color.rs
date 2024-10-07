use rand::prelude::*;
use rand_distr::{Distribution, Normal};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct RGBColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl RGBColor {
    pub fn new() -> RGBColor {
        RGBColor { r: 0, g: 0, b: 0 }
    }

    pub fn from_normalized(r: f32, g: f32, b: f32) -> RGBColor {
        RGBColor {
            r: (r * 255.).round() as u8,
            g: (g * 255.).round() as u8,
            b: (b * 255.).round() as u8,
        }
    }

    pub fn normalize(&self) -> SGRBColor {
        SGRBColor {
            r: self.r as f32 / 255.,
            g: self.g as f32 / 255.,
            b: self.b as f32 / 255.,
        }
    }

    pub fn random() -> RGBColor {
        let color = RGBColor {
            r: random::<u8>(),
            g: random::<u8>(),
            b: random::<u8>(),
        };
        color
    }

    pub fn from_linear(linear_sgrb: SGRBColor) -> RGBColor {
        let mut values: [f32; 3] = [linear_sgrb.r, linear_sgrb.g, linear_sgrb.b];

        for i in 0..3 {
            let mut val = values[i];
            if val >= 0.0031308 {
                val = 1.055 * val.powf(1.0 / 2.4) - 0.055;
            } else {
                val = 12.92 * val;
            }
            values[i] = val;
        }

        RGBColor::from_normalized(values[0], values[1], values[2])
    }

    pub fn to_linear(&self) -> SGRBColor {
        let normalized = self.normalize();
        let mut values: [f32; 3] = [normalized.r, normalized.g, normalized.b];

        for i in 0..3 {
            let mut val = values[i];
            if val >= 0.04045 {
                val = ((val + 0.055) / (1. + 0.055)).powf(2.4);
            } else {
                val = val / 12.92;
            }
            values[i] = val;
        }

        SGRBColor {
            r: values[0],
            g: values[1],
            b: values[2],
        }
    }

    pub fn from_oklab(oklab: OklabColor) -> RGBColor {
        RGBColor::from_linear(SGRBColor::from_oklab(oklab))
    }

    pub fn to_oklab(&self) -> OklabColor {
        self.to_linear().to_oklab()
    }
}

#[derive(Copy, Clone, Debug)]
pub struct SGRBColor {
    pub r: f32,
    pub g: f32,
    pub b: f32,
}

impl SGRBColor {
    pub fn new() -> SGRBColor {
        SGRBColor {
            r: 0.,
            g: 0.,
            b: 0.,
        }
    }

    // Conversion codes are adapted from Björn Ottosson's work. (/licenses/oklab)
    pub fn from_oklab(oklab: OklabColor) -> SGRBColor {
        let l = oklab.L + 0.3963377774 * oklab.a + 0.2158037573 * oklab.b;
        let m = oklab.L - 0.1055613458 * oklab.a - 0.0638541728 * oklab.b;
        let s = oklab.L - 0.0894841775 * oklab.a - 1.2914855480 * oklab.b;

        let l = l * l * l;
        let m = m * m * m;
        let s = s * s * s;

        let r = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s;
        let g = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s;
        let b = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s;

        SGRBColor { r, g, b }
    }

    pub fn to_oklab(&self) -> OklabColor {
        let l = 0.4122214708 * self.r + 0.5363325363 * self.g + 0.0514459929 * self.b;
        let m = 0.2119034982 * self.r + 0.6806995451 * self.g + 0.1073969566 * self.b;
        let s = 0.0883024619 * self.r + 0.2817188376 * self.g + 0.6299787005 * self.b;

        let l = l.cbrt();
        let m = m.cbrt();
        let s = s.cbrt();

        OklabColor {
            L: 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            a: 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            b: 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
        }
    }

    // Clip RGB values to lie between [0..1]
    pub fn clip(&mut self) {
        self.r = self.r.clamp(0.0, 1.0);
        self.g = self.g.clamp(0.0, 1.0);
        self.b = self.b.clamp(0.0, 1.0);
    }

    pub fn is_valid(&self) -> bool {
        let values = [self.r, self.g, self.b];
        for val in values {
            if val < 0. || val > 1. {
                return false;
            }
        }
        return true;
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
#[allow(non_snake_case)]
pub struct OklabColor {
    pub L: f32,
    pub a: f32,
    pub b: f32,
}

pub enum OklabResult {
    Ok(OklabColor),
    Clipped(OklabColor), // returned if color was clipped to be in SGRB
}

impl OklabResult {
    pub fn unwrap(self) -> OklabColor {
        match self {
            Self::Ok(c) => c,
            Self::Clipped(c) => c,
        }
    }
}

impl OklabColor {
    pub fn new() -> OklabColor {
        OklabColor {
            L: 0.,
            a: 0.,
            b: 0.,
        }
    }

    #[allow(non_snake_case)]
    pub fn random(
        min_L: f32,
        max_L: f32,
        min_a: f32,
        max_a: f32,
        min_b: f32,
        max_b: f32,
    ) -> OklabColor {
        OklabColor {
            L: thread_rng().gen_range(min_L..max_L),
            a: thread_rng().gen_range(min_a..max_a),
            b: thread_rng().gen_range(min_b..max_b),
        }
    }

    // Means: [ 0.63743672  0.00487444 -0.00147413]
    // Stds: [0.16545165 0.11717777 0.11457232]
    // maxs: [0.99999999 0.2762164  0.19856975]
    // mins: [ 0.         -0.23388757 -0.31152815]
    #[allow(non_snake_case)]
    pub fn random_in_sgrb(min_L: f32, max_L: f32, max_tries: u32) -> OklabResult {
        const MIN_A: f32 = -0.23388757;
        const MAX_A: f32 = 0.2762164;
        const MIN_B: f32 = -0.31152815;
        const MAX_B: f32 = 0.19856975;

        let mut tries = 0;
        let mut color = Self::random(min_L, max_L, MIN_A, MAX_A, MIN_B, MAX_B);
        let mut fitness = color.sgrb_fitness();
        let mut best_color = color;
        let mut best_fitness = fitness;

        while fitness < 0f32 {
            if tries == max_tries {
                let mut clipped_sgrb = SGRBColor::from_oklab(best_color);
                clipped_sgrb.clip();
                return OklabResult::Clipped(clipped_sgrb.to_oklab());
            }
            color = Self::random(min_L, max_L, MIN_A, MAX_A, MIN_B, MAX_B);
            fitness = color.sgrb_fitness();
            if fitness > best_fitness {
                best_fitness = fitness;
                best_color = color;
            }
            tries += 1;
        }

        return OklabResult::Ok(color);
    }

    #[allow(non_snake_case)]
    fn step(&self, delta_L: f32, delta_a: f32, delta_b: f32) -> OklabColor {
        OklabColor {
            L: (self.L + delta_L).clamp(0.0, 1.0),
            a: self.a + delta_a,
            b: self.b + delta_b,
        }
    }

    fn random_step(&self, step_size: f32) -> OklabColor {
        let normal_dist = Normal::new(0.0, 1.0).expect("Failed to create a normal distribution.");

        let mut step = [0f32; 3];
        let mut norm = 0.0;
        for x in &mut step {
            *x = normal_dist.sample(&mut rand::thread_rng());
            norm += (*x).powf(2.0);
        }
        norm = norm.sqrt();
        for x in &mut step {
            *x = *x / norm * step_size;
        }

        self.step(step[0], step[1], step[2])
    }

    // take a random step and return a color that is valid in SGRB space
    // If this fails after max_tries, it clips and returns the closest color
    pub fn random_step_sgrb(&self, step_size: f32, max_tries: u32) -> OklabResult {
        let mut best_color: OklabColor = self.random_step(step_size);
        let mut best_fitness = best_color.sgrb_fitness();
        if best_fitness >= 0.0 {
            return OklabResult::Ok(best_color);
        }

        let mut num_tries = 1u32;
        while num_tries < max_tries {
            let color = self.random_step(step_size);
            let fitness = color.sgrb_fitness();

            if fitness >= 0.0 {
                return OklabResult::Ok(color);
            } else if fitness > best_fitness {
                best_fitness = fitness;
                best_color = color;
            }
            num_tries += 1;
        }

        let mut clipped_sgrb = SGRBColor::from_oklab(best_color);
        clipped_sgrb.clip();
        OklabResult::Clipped(clipped_sgrb.to_oklab())
    }

    // Return the euclidian distance between two Oklab colors
    pub fn distance(&self, other: &OklabColor) -> f32 {
        let dist = (self.L - other.L).powf(2.0)
            + (self.a - other.a).powf(2.0)
            + (self.b - other.b).powf(2.0);

        dist.sqrt()
    }

    // Return the squared distance between two Oklab colors
    pub fn squared_distance(&self, other: &OklabColor) -> f32 {
        let dist_sq = (self.L - other.L).powf(2.0)
            + (self.a - other.a).powf(2.0)
            + (self.b - other.b).powf(2.0);

        dist_sq
    }

    // Returns a value signifying closeness to a valid SGRB color.
    // If the color can be converted to SGRB, this function returns 0.
    // Otherwise, it returns a negative value.
    fn sgrb_fitness(&self) -> f32 {
        let mut fitness = 0f32;
        let sgrb = SGRBColor::from_oklab(*self);

        let values: [f32; 3] = [sgrb.r, sgrb.g, sgrb.b];

        for val in &values {
            if *val > 1.0 {
                fitness -= *val - 1.0;
            } else if *val < 0.0 {
                fitness += *val;
            }
        }

        return fitness;
    }
}

#[cfg(test)]
mod tests {
    use float_cmp::approx_eq;

    use super::*;

    #[test]
    fn linear_to_rgb() {
        let converted = RGBColor::from_linear(SGRBColor {
            r: (0.5),
            g: (0.3),
            b: (0.87),
        });
        assert_eq!(
            converted,
            RGBColor {
                r: 188,
                g: 149,
                b: 240
            }
        );

        let converted = RGBColor::from_linear(SGRBColor {
            r: (0.),
            g: (1.),
            b: (0.25),
        });
        assert_eq!(
            converted,
            RGBColor {
                r: 0,
                g: 255,
                b: 137
            }
        );
    }

    #[test]
    fn rgb_to_linear() {
        let rgb = RGBColor {
            r: 50,
            g: 149,
            b: 255,
        };
        let linear = rgb.to_linear();

        assert!(approx_eq!(f32, 0.031896033, linear.r, ulps = 2));
        assert!(approx_eq!(f32, 0.30054379, linear.g, ulps = 2));
        assert!(approx_eq!(f32, 1., linear.b, ulps = 2));

        let rgb = RGBColor { r: 0, g: 8, b: 100 };
        let linear = rgb.to_linear();

        assert!(approx_eq!(f32, 0., linear.r, ulps = 2));
        assert!(approx_eq!(f32, 0.0024282159, linear.g, ulps = 2));
        assert!(approx_eq!(f32, 0.12743768, linear.b, ulps = 2));
    }

    #[test]
    fn rgb_to_oklab() {
        let rgb = RGBColor {
            r: 166,
            g: 102,
            b: 40,
        };
        let oklab = rgb.to_oklab();

        // reference values taken from css oklab().
        // The values don't match closely, probably due to Björn Ottosson's updated matrix values in 2021.
        assert!(approx_eq!(f32, 0.57, oklab.L, epsilon = 0.001));
        assert!(approx_eq!(f32, 0.055, oklab.a, epsilon = 0.001));
        assert!(approx_eq!(f32, 0.098, oklab.b, epsilon = 0.001));
    }

    #[test]
    fn oklab_to_rgb() {
        let oklab = OklabColor {
            L: 0.57,
            a: 0.055,
            b: 0.098,
        };
        let rgb = RGBColor::from_oklab(oklab);

        dbg!(oklab);
        dbg!(rgb);
        assert_eq!(
            rgb,
            RGBColor {
                r: 166,
                g: 102,
                b: 40
            }
        );
    }
}
