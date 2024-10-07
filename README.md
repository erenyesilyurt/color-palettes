# Color Palettes

Generate colors with a perceptually uniform distance between them.

```rust
use color_palettes::*;

fn main() {
    let num_colors = 16;
    let config = GeneratorConfig::default();
    let palette = generate_palette(num_colors, &config).unwrap();
    save_palette(&palette, "palette.png");
}
```

You can adjust the distance between subsequent colors in the palette by changing similarity:

```rust
    let mut config = GeneratorConfig::default();
    config.similarity = 0.8; // takes values between 0 and 1
    let palette = generate_palette(num_colors, &config).unwrap();
```