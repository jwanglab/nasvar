//! PNG rendering backend with built-in rasterization and encoding.

use super::Backend;
use crate::plotting::style::{DominantBaseline, FillStyle, LineStyle, Marker, MarkerStyle, TextAnchor, TextStyle};

// ─── Bitmap font (5x7, ASCII 32–126) ────────────────────────────────────────

/// 5-pixel-wide, 7-pixel-tall bitmap font covering ASCII 32–126.
/// Each character is stored as 7 bytes; bits 4..0 of each byte are the 5 pixel columns
/// (bit 4 = leftmost column).
const FONT: [[u8; 7]; 95] = [
    [0x00,0x00,0x00,0x00,0x00,0x00,0x00], // 32 ' '
    [0x04,0x04,0x04,0x04,0x04,0x00,0x04], // 33 '!'
    [0x0A,0x0A,0x00,0x00,0x00,0x00,0x00], // 34 '"'
    [0x0A,0x1F,0x0A,0x0A,0x1F,0x0A,0x00], // 35 '#'
    [0x04,0x0F,0x14,0x0E,0x05,0x1E,0x04], // 36 '$'
    [0x19,0x19,0x02,0x04,0x08,0x13,0x13], // 37 '%'
    [0x0C,0x12,0x0C,0x0D,0x12,0x12,0x0D], // 38 '&'
    [0x04,0x04,0x00,0x00,0x00,0x00,0x00], // 39 '''
    [0x02,0x04,0x08,0x08,0x08,0x04,0x02], // 40 '('
    [0x08,0x04,0x02,0x02,0x02,0x04,0x08], // 41 ')'
    [0x00,0x04,0x15,0x0E,0x15,0x04,0x00], // 42 '*'
    [0x00,0x04,0x04,0x1F,0x04,0x04,0x00], // 43 '+'
    [0x00,0x00,0x00,0x00,0x00,0x04,0x08], // 44 ','
    [0x00,0x00,0x00,0x1F,0x00,0x00,0x00], // 45 '-'
    [0x00,0x00,0x00,0x00,0x00,0x00,0x04], // 46 '.'
    [0x01,0x01,0x02,0x04,0x08,0x10,0x10], // 47 '/'
    [0x0E,0x11,0x13,0x15,0x19,0x11,0x0E], // 48 '0'
    [0x04,0x0C,0x04,0x04,0x04,0x04,0x0E], // 49 '1'
    [0x0E,0x11,0x01,0x06,0x08,0x10,0x1F], // 50 '2'
    [0x0E,0x11,0x01,0x06,0x01,0x11,0x0E], // 51 '3'
    [0x02,0x06,0x0A,0x12,0x1F,0x02,0x02], // 52 '4'
    [0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E], // 53 '5'
    [0x06,0x08,0x10,0x1E,0x11,0x11,0x0E], // 54 '6'
    [0x1F,0x01,0x02,0x04,0x08,0x08,0x08], // 55 '7'
    [0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E], // 56 '8'
    [0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C], // 57 '9'
    [0x00,0x00,0x04,0x00,0x00,0x04,0x00], // 58 ':'
    [0x00,0x00,0x04,0x00,0x00,0x04,0x08], // 59 ';'
    [0x02,0x04,0x08,0x10,0x08,0x04,0x02], // 60 '<'
    [0x00,0x00,0x1F,0x00,0x1F,0x00,0x00], // 61 '='
    [0x08,0x04,0x02,0x01,0x02,0x04,0x08], // 62 '>'
    [0x0E,0x11,0x01,0x02,0x04,0x00,0x04], // 63 '?'
    [0x0E,0x11,0x17,0x15,0x17,0x10,0x0E], // 64 '@'
    [0x0E,0x11,0x11,0x1F,0x11,0x11,0x11], // 65 'A'
    [0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E], // 66 'B'
    [0x0E,0x11,0x10,0x10,0x10,0x11,0x0E], // 67 'C'
    [0x1E,0x11,0x11,0x11,0x11,0x11,0x1E], // 68 'D'
    [0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F], // 69 'E'
    [0x1F,0x10,0x10,0x1E,0x10,0x10,0x10], // 70 'F'
    [0x0E,0x11,0x10,0x17,0x11,0x11,0x0F], // 71 'G'
    [0x11,0x11,0x11,0x1F,0x11,0x11,0x11], // 72 'H'
    [0x0E,0x04,0x04,0x04,0x04,0x04,0x0E], // 73 'I'
    [0x07,0x02,0x02,0x02,0x02,0x12,0x0C], // 74 'J'
    [0x11,0x12,0x14,0x18,0x14,0x12,0x11], // 75 'K'
    [0x10,0x10,0x10,0x10,0x10,0x10,0x1F], // 76 'L'
    [0x11,0x1B,0x15,0x15,0x11,0x11,0x11], // 77 'M'
    [0x11,0x19,0x15,0x13,0x11,0x11,0x11], // 78 'N'
    [0x0E,0x11,0x11,0x11,0x11,0x11,0x0E], // 79 'O'
    [0x1E,0x11,0x11,0x1E,0x10,0x10,0x10], // 80 'P'
    [0x0E,0x11,0x11,0x11,0x15,0x12,0x0D], // 81 'Q'
    [0x1E,0x11,0x11,0x1E,0x14,0x12,0x11], // 82 'R'
    [0x0E,0x11,0x10,0x0E,0x01,0x11,0x0E], // 83 'S'
    [0x1F,0x04,0x04,0x04,0x04,0x04,0x04], // 84 'T'
    [0x11,0x11,0x11,0x11,0x11,0x11,0x0E], // 85 'U'
    [0x11,0x11,0x11,0x11,0x0A,0x0A,0x04], // 86 'V'
    [0x11,0x11,0x11,0x15,0x15,0x1B,0x11], // 87 'W'
    [0x11,0x11,0x0A,0x04,0x0A,0x11,0x11], // 88 'X'
    [0x11,0x11,0x0A,0x04,0x04,0x04,0x04], // 89 'Y'
    [0x1F,0x01,0x02,0x04,0x08,0x10,0x1F], // 90 'Z'
    [0x0E,0x08,0x08,0x08,0x08,0x08,0x0E], // 91 '['
    [0x10,0x10,0x08,0x04,0x02,0x01,0x01], // 92 '\'
    [0x0E,0x02,0x02,0x02,0x02,0x02,0x0E], // 93 ']'
    [0x04,0x0A,0x11,0x00,0x00,0x00,0x00], // 94 '^'
    [0x00,0x00,0x00,0x00,0x00,0x00,0x1F], // 95 '_'
    [0x08,0x04,0x00,0x00,0x00,0x00,0x00], // 96 '`'
    [0x00,0x00,0x0E,0x01,0x0F,0x11,0x0F], // 97 'a'
    [0x10,0x10,0x1E,0x11,0x11,0x11,0x1E], // 98 'b'
    [0x00,0x00,0x0E,0x11,0x10,0x11,0x0E], // 99 'c'
    [0x01,0x01,0x0F,0x11,0x11,0x11,0x0F], // 100 'd'
    [0x00,0x00,0x0E,0x11,0x1F,0x10,0x0E], // 101 'e'
    [0x06,0x08,0x1E,0x08,0x08,0x08,0x08], // 102 'f'
    [0x00,0x00,0x0F,0x11,0x0F,0x01,0x0E], // 103 'g'
    [0x10,0x10,0x1E,0x11,0x11,0x11,0x11], // 104 'h'
    [0x04,0x00,0x0C,0x04,0x04,0x04,0x0E], // 105 'i'
    [0x02,0x00,0x06,0x02,0x02,0x12,0x0C], // 106 'j'
    [0x10,0x10,0x12,0x14,0x18,0x14,0x12], // 107 'k'
    [0x0C,0x04,0x04,0x04,0x04,0x04,0x0E], // 108 'l'
    [0x00,0x00,0x1A,0x15,0x15,0x15,0x15], // 109 'm'
    [0x00,0x00,0x1E,0x11,0x11,0x11,0x11], // 110 'n'
    [0x00,0x00,0x0E,0x11,0x11,0x11,0x0E], // 111 'o'
    [0x00,0x00,0x1E,0x11,0x1E,0x10,0x10], // 112 'p'
    [0x00,0x00,0x0F,0x11,0x0F,0x01,0x01], // 113 'q'
    [0x00,0x00,0x16,0x19,0x10,0x10,0x10], // 114 'r'
    [0x00,0x00,0x0F,0x10,0x0E,0x01,0x1E], // 115 's'
    [0x08,0x08,0x1E,0x08,0x08,0x09,0x06], // 116 't'
    [0x00,0x00,0x11,0x11,0x11,0x13,0x0D], // 117 'u'
    [0x00,0x00,0x11,0x11,0x11,0x0A,0x04], // 118 'v'
    [0x00,0x00,0x11,0x11,0x15,0x15,0x0A], // 119 'w'
    [0x00,0x00,0x11,0x0A,0x04,0x0A,0x11], // 120 'x'
    [0x00,0x00,0x11,0x11,0x0F,0x01,0x0E], // 121 'y'
    [0x00,0x00,0x1F,0x02,0x04,0x08,0x1F], // 122 'z'
    [0x02,0x04,0x04,0x08,0x04,0x04,0x02], // 123 '{'
    [0x04,0x04,0x04,0x04,0x04,0x04,0x04], // 124 '|'
    [0x08,0x04,0x04,0x02,0x04,0x04,0x08], // 125 '}'
    [0x00,0x0D,0x12,0x00,0x00,0x00,0x00], // 126 '~'
];

const FONT_W: usize = 5;
const FONT_H: usize = 7;
const FONT_SPACING: usize = 1; // extra pixels between characters

// ─── Clipping ────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy)]
struct ClipRect {
    x: f64,
    y: f64,
    w: f64,
    h: f64,
}

// ─── PngBackend ──────────────────────────────────────────────────────────────

/// A raster rendering backend that produces PNG output.
///
/// Implements all drawing primitives with anti-aliased rendering and
/// includes a minimal built-in PNG encoder (no external dependencies).
pub struct PngBackend {
    /// Image width in pixels
    pub width: u32,
    /// Image height in pixels
    pub height: u32,
    /// RGBA pixel buffer (4 bytes per pixel, row-major)
    pixels: Vec<u8>,
    /// Stack of active clip rectangles
    clip_stack: Vec<ClipRect>,
}

impl PngBackend {
    /// Create a new PNG backend with the given dimensions.
    /// Initializes with a white background.
    pub fn new(width: u32, height: u32) -> Self {
        // Initialize to white background (RGBA: 255, 255, 255, 255)
        let pixels = vec![255u8; (width * height * 4) as usize];
        PngBackend {
            width,
            height,
            pixels,
            clip_stack: Vec::new(),
        }
    }

    /// Encode the pixel buffer as a PNG file.
    pub fn encode_png(&self) -> Vec<u8> {
        encode_png_data(self.width, self.height, &self.pixels)
    }

    // ── Pixel-level helpers ──────────────────────────────────────────────

    fn in_bounds(&self, x: i32, y: i32) -> bool {
        x >= 0 && y >= 0 && (x as u32) < self.width && (y as u32) < self.height
    }

    fn in_clip(&self, x: f64, y: f64) -> bool {
        if let Some(c) = self.clip_stack.last() {
            x >= c.x && x < c.x + c.w && y >= c.y && y < c.y + c.h
        } else {
            true
        }
    }

    /// Blend a single pixel with alpha compositing (source-over).
    fn blend_pixel(&mut self, x: i32, y: i32, r: u8, g: u8, b: u8, a: u8) {
        if !self.in_bounds(x, y) || !self.in_clip(x as f64, y as f64) {
            return;
        }
        if a == 0 {
            return;
        }
        let idx = ((y as u32 * self.width + x as u32) * 4) as usize;
        if a == 255 {
            self.pixels[idx] = r;
            self.pixels[idx + 1] = g;
            self.pixels[idx + 2] = b;
            self.pixels[idx + 3] = 255;
            return;
        }
        let sa = a as f64 / 255.0;
        let da = self.pixels[idx + 3] as f64 / 255.0;
        let out_a = sa + da * (1.0 - sa);
        if out_a <= 0.0 {
            return;
        }
        let blend = |s: u8, d: u8| -> u8 {
            ((s as f64 * sa + d as f64 * da * (1.0 - sa)) / out_a).round() as u8
        };
        self.pixels[idx] = blend(r, self.pixels[idx]);
        self.pixels[idx + 1] = blend(g, self.pixels[idx + 1]);
        self.pixels[idx + 2] = blend(b, self.pixels[idx + 2]);
        self.pixels[idx + 3] = (out_a * 255.0).round() as u8;
    }

    // ── Shape primitives ─────────────────────────────────────────────────

    fn fill_rect_px(&mut self, x: f64, y: f64, w: f64, h: f64, r: u8, g: u8, b: u8, a: u8) {
        let x0 = x.floor() as i32;
        let y0 = y.floor() as i32;
        let x1 = (x + w).ceil() as i32;
        let y1 = (y + h).ceil() as i32;
        for py in y0..y1 {
            for px in x0..x1 {
                self.blend_pixel(px, py, r, g, b, a);
            }
        }
    }

    fn stroke_rect_px(&mut self, x: f64, y: f64, w: f64, h: f64, r: u8, g: u8, b: u8, a: u8, sw: f64) {
        // top
        self.fill_rect_px(x, y, w, sw, r, g, b, a);
        // bottom
        self.fill_rect_px(x, y + h - sw, w, sw, r, g, b, a);
        // left
        self.fill_rect_px(x, y, sw, h, r, g, b, a);
        // right
        self.fill_rect_px(x + w - sw, y, sw, h, r, g, b, a);
    }

    fn fill_circle_px(&mut self, cx: f64, cy: f64, radius: f64, r: u8, g: u8, b: u8, a: u8) {
        let r2 = radius * radius;
        let x0 = (cx - radius - 1.0).floor() as i32;
        let y0 = (cy - radius - 1.0).floor() as i32;
        let x1 = (cx + radius + 1.0).ceil() as i32;
        let y1 = (cy + radius + 1.0).ceil() as i32;
        for py in y0..=y1 {
            for px in x0..=x1 {
                let dx = px as f64 + 0.5 - cx;
                let dy = py as f64 + 0.5 - cy;
                let dist2 = dx * dx + dy * dy;
                if dist2 <= r2 {
                    self.blend_pixel(px, py, r, g, b, a);
                } else if dist2 <= (radius + 1.0) * (radius + 1.0) {
                    // Anti-alias at edge
                    let dist = dist2.sqrt();
                    let coverage = (radius + 0.5 - dist).clamp(0.0, 1.0);
                    let aa = (a as f64 * coverage).round() as u8;
                    if aa > 0 {
                        self.blend_pixel(px, py, r, g, b, aa);
                    }
                }
            }
        }
    }

    fn stroke_circle_px(&mut self, cx: f64, cy: f64, radius: f64, r: u8, g: u8, b: u8, a: u8, sw: f64) {
        let outer = radius + sw / 2.0;
        let inner = (radius - sw / 2.0).max(0.0);
        let x0 = (cx - outer - 1.0).floor() as i32;
        let y0 = (cy - outer - 1.0).floor() as i32;
        let x1 = (cx + outer + 1.0).ceil() as i32;
        let y1 = (cy + outer + 1.0).ceil() as i32;
        for py in y0..=y1 {
            for px in x0..=x1 {
                let dx = px as f64 + 0.5 - cx;
                let dy = py as f64 + 0.5 - cy;
                let dist = (dx * dx + dy * dy).sqrt();
                if dist >= inner - 0.5 && dist <= outer + 0.5 {
                    let coverage_outer = (outer + 0.5 - dist).clamp(0.0, 1.0);
                    let coverage_inner = (dist - inner + 0.5).clamp(0.0, 1.0);
                    let coverage = coverage_outer.min(coverage_inner);
                    let aa = (a as f64 * coverage).round() as u8;
                    if aa > 0 {
                        self.blend_pixel(px, py, r, g, b, aa);
                    }
                }
            }
        }
    }

    /// Draw an anti-aliased line using Xiaolin Wu's algorithm.
    fn draw_line_aa(&mut self, mut x0: f64, mut y0: f64, mut x1: f64, mut y1: f64,
                     r: u8, g: u8, b: u8, a: u8) {
        let steep = (y1 - y0).abs() > (x1 - x0).abs();
        if steep {
            std::mem::swap(&mut x0, &mut y0);
            std::mem::swap(&mut x1, &mut y1);
        }
        if x0 > x1 {
            std::mem::swap(&mut x0, &mut x1);
            std::mem::swap(&mut y0, &mut y1);
        }

        let dx = x1 - x0;
        let dy = y1 - y0;
        let gradient = if dx.abs() < 1e-10 { 1.0 } else { dy / dx };

        // Helper to plot with intensity
        let plot = |backend: &mut Self, x: i32, y: i32, intensity: f64| {
            let aa = (a as f64 * intensity).round() as u8;
            if aa > 0 {
                if steep {
                    backend.blend_pixel(y, x, r, g, b, aa);
                } else {
                    backend.blend_pixel(x, y, r, g, b, aa);
                }
            }
        };

        // First endpoint
        let xend = x0.round();
        let yend = y0 + gradient * (xend - x0);
        let xgap = 1.0 - (x0 + 0.5).fract();
        let xpxl1 = xend as i32;
        let ypxl1 = yend.floor() as i32;
        plot(self, xpxl1, ypxl1, (1.0 - yend.fract()) * xgap);
        plot(self, xpxl1, ypxl1 + 1, yend.fract() * xgap);
        let mut intery = yend + gradient;

        // Second endpoint
        let xend = x1.round();
        let yend = y1 + gradient * (xend - x1);
        let xgap = (x1 + 0.5).fract();
        let xpxl2 = xend as i32;
        let ypxl2 = yend.floor() as i32;
        plot(self, xpxl2, ypxl2, (1.0 - yend.fract()) * xgap);
        plot(self, xpxl2, ypxl2 + 1, yend.fract() * xgap);

        // Main loop
        for x in (xpxl1 + 1)..xpxl2 {
            let y_floor = intery.floor() as i32;
            let frac = intery.fract();
            plot(self, x, y_floor, 1.0 - frac);
            plot(self, x, y_floor + 1, frac);
            intery += gradient;
        }
    }

    /// Draw a thick line with anti-aliased edges.
    fn draw_thick_line(&mut self, x0: f64, y0: f64, x1: f64, y1: f64, width: f64,
                        r: u8, g: u8, b: u8, a: u8) {
        // For thin lines, use anti-aliased line drawing
        if width <= 1.0 {
            self.draw_line_aa(x0, y0, x1, y1, r, g, b, a);
            return;
        }

        let dx = x1 - x0;
        let dy = y1 - y0;
        let len = (dx * dx + dy * dy).sqrt();
        if len < 1e-10 {
            // Just a point - draw a small filled circle
            self.fill_circle_px(x0, y0, width / 2.0, r, g, b, a);
            return;
        }

        // For thick lines, use distance-based anti-aliasing
        let hw = width / 2.0;

        // Line direction
        let lx = dx / len;
        let ly = dy / len;

        // Bounding box with padding for AA
        let padding = hw + 1.5;
        let min_x = (x0.min(x1) - padding).floor() as i32;
        let max_x = (x0.max(x1) + padding).ceil() as i32;
        let min_y = (y0.min(y1) - padding).floor() as i32;
        let max_y = (y0.max(y1) + padding).ceil() as i32;

        for py_i in min_y..=max_y {
            for px_i in min_x..=max_x {
                let px_f = px_i as f64 + 0.5;
                let py_f = py_i as f64 + 0.5;

                // Vector from line start to pixel
                let vx = px_f - x0;
                let vy = py_f - y0;

                // Project onto line direction to get parameter t
                let t = vx * lx + vy * ly;

                // Clamp t to line segment
                let t_clamped = t.clamp(0.0, len);

                // Closest point on line segment
                let closest_x = x0 + lx * t_clamped;
                let closest_y = y0 + ly * t_clamped;

                // Distance from pixel to closest point
                let dist_x = px_f - closest_x;
                let dist_y = py_f - closest_y;
                let dist = (dist_x * dist_x + dist_y * dist_y).sqrt();

                // Anti-aliased coverage based on distance
                if dist < hw + 1.0 {
                    let coverage = if dist <= hw - 0.5 {
                        1.0
                    } else {
                        (hw + 0.5 - dist).clamp(0.0, 1.0)
                    };
                    let aa = (a as f64 * coverage).round() as u8;
                    if aa > 0 {
                        self.blend_pixel(px_i, py_i, r, g, b, aa);
                    }
                }
            }
        }
    }

    /// Draw a line with a dash pattern.
    fn draw_dashed_line(&mut self, x0: f64, y0: f64, x1: f64, y1: f64,
                         width: f64, pattern: &[f64], r: u8, g: u8, b: u8, a: u8) {
        let dx = x1 - x0;
        let dy = y1 - y0;
        let total_len = (dx * dx + dy * dy).sqrt();
        if total_len < 1e-10 || pattern.is_empty() {
            return;
        }
        let ux = dx / total_len;
        let uy = dy / total_len;
        let mut dist = 0.0;
        let mut idx = 0;
        let mut drawing = true;
        while dist < total_len {
            let seg_len = pattern[idx % pattern.len()].min(total_len - dist);
            if drawing {
                let sx = x0 + ux * dist;
                let sy = y0 + uy * dist;
                let ex = x0 + ux * (dist + seg_len);
                let ey = y0 + uy * (dist + seg_len);
                self.draw_thick_line(sx, sy, ex, ey, width, r, g, b, a);
            }
            dist += seg_len;
            idx += 1;
            drawing = !drawing;
        }
    }

    /// Fill a polygon using scanline algorithm.
    fn fill_polygon_px(&mut self, points: &[(f64, f64)], r: u8, g: u8, b: u8, a: u8) {
        if points.len() < 3 {
            return;
        }
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        for &(_, y) in points {
            min_y = min_y.min(y);
            max_y = max_y.max(y);
        }
        let y_start = min_y.floor() as i32;
        let y_end = max_y.ceil() as i32;
        let n = points.len();

        for y in y_start..=y_end {
            let yf = y as f64 + 0.5;
            let mut intersections = Vec::new();
            for i in 0..n {
                let j = (i + 1) % n;
                let (_, y0) = points[i];
                let (_, y1) = points[j];
                if (y0 <= yf && y1 > yf) || (y1 <= yf && y0 > yf) {
                    let (x0, _) = points[i];
                    let (x1, _) = points[j];
                    let t = (yf - y0) / (y1 - y0);
                    let x = x0 + t * (x1 - x0);
                    intersections.push(x);
                }
            }
            intersections.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for pair in intersections.chunks(2) {
                if pair.len() == 2 {
                    let x_start = pair[0].floor() as i32;
                    let x_end = pair[1].ceil() as i32;
                    for x in x_start..=x_end {
                        self.blend_pixel(x, y, r, g, b, a);
                    }
                }
            }
        }
    }

    // ── Text rendering ───────────────────────────────────────────────────

    /// Measure text width at the given font scale.
    fn text_width(&self, text: &str, scale: usize) -> f64 {
        let char_w = FONT_W * scale;
        let spacing = FONT_SPACING * scale;
        let n = text.chars().filter(|c| *c >= ' ' && *c <= '~').count();
        if n == 0 {
            return 0.0;
        }
        (n * char_w + (n - 1) * spacing) as f64
    }

    /// Render a text string. Rotation in degrees.
    fn render_text(&mut self, x: f64, y: f64, text: &str, scale: usize,
                    anchor: &TextAnchor, baseline: &DominantBaseline,
                    rotation: f64, r: u8, g: u8, b: u8, a: u8) {
        let total_w = self.text_width(text, scale);
        let total_h = (FONT_H * scale) as f64;

        // Horizontal offset based on anchor
        let x_offset = match anchor {
            TextAnchor::Start => 0.0,
            TextAnchor::Middle => -total_w / 2.0,
            TextAnchor::End => -total_w,
        };

        // Vertical offset based on baseline
        let y_offset = match baseline {
            DominantBaseline::Auto | DominantBaseline::Alphabetic => -total_h,
            DominantBaseline::Middle => -total_h / 2.0,
            DominantBaseline::Hanging => 0.0,
            DominantBaseline::Ideographic => -total_h,
        };

        let use_rotation = rotation.abs() > 0.5;
        let cos_a = if use_rotation { rotation.to_radians().cos() } else { 1.0 };
        let sin_a = if use_rotation { rotation.to_radians().sin() } else { 0.0 };

        let char_w = FONT_W * scale;
        let spacing = FONT_SPACING * scale;
        let mut cx = 0.0;

        for ch in text.chars() {
            let idx = ch as usize;
            if idx < 32 || idx > 126 {
                cx += (char_w + spacing) as f64;
                continue;
            }
            let glyph = &FONT[idx - 32];
            for row in 0..FONT_H {
                let bits = glyph[row];
                for col in 0..FONT_W {
                    if bits & (1 << (4 - col)) != 0 {
                        // For each set bit, fill a scale×scale block
                        for sy in 0..scale {
                            for sx in 0..scale {
                                let lx = cx + (col * scale + sx) as f64 + x_offset;
                                let ly = (row * scale + sy) as f64 + y_offset;
                                let (px, py) = if use_rotation {
                                    (x + lx * cos_a - ly * sin_a,
                                     y + lx * sin_a + ly * cos_a)
                                } else {
                                    (x + lx, y + ly)
                                };
                                self.blend_pixel(px.round() as i32, py.round() as i32, r, g, b, a);
                            }
                        }
                    }
                }
            }
            cx += (char_w + spacing) as f64;
        }
    }

    // ── Marker rendering ─────────────────────────────────────────────────

    fn render_marker_at(&mut self, cx: f64, cy: f64, style: &MarkerStyle) {
        let (fr, fg, fb, fa) = style.fill.to_rgba();
        let fa = (fa as f64 * style.fill_opacity).round() as u8;
        let (er, eg, eb, ea) = style.edge_color.to_rgba();
        let radius = style.size / 2.0;

        match &style.marker {
            Marker::None => {}
            Marker::Circle => {
                self.fill_circle_px(cx, cy, radius, fr, fg, fb, fa);
                if style.edge_width > 0.0 {
                    self.stroke_circle_px(cx, cy, radius, er, eg, eb, ea, style.edge_width);
                }
            }
            Marker::Cross => {
                let s = radius * 0.707;
                let line_w = style.edge_width.max(1.5);
                self.draw_thick_line(cx - s, cy - s, cx + s, cy + s, line_w, er, eg, eb, ea);
                self.draw_thick_line(cx - s, cy + s, cx + s, cy - s, line_w, er, eg, eb, ea);
            }
            marker => {
                // Build polygon points from the marker type
                if let Some(points) = marker_polygon_points(marker, radius) {
                    let translated: Vec<(f64, f64)> = points.iter()
                        .map(|(x, y)| (cx + x, cy + y))
                        .collect();
                    self.fill_polygon_px(&translated, fr, fg, fb, fa);
                    if style.edge_width > 0.0 {
                        let n = translated.len();
                        for i in 0..n {
                            let j = (i + 1) % n;
                            self.draw_thick_line(
                                translated[i].0, translated[i].1,
                                translated[j].0, translated[j].1,
                                style.edge_width, er, eg, eb, ea,
                            );
                        }
                    }
                }
            }
        }
    }
}

/// Compute polygon vertices for a marker shape centered at origin.
fn marker_polygon_points(marker: &Marker, size: f64) -> Option<Vec<(f64, f64)>> {
    match marker {
        Marker::None | Marker::Circle | Marker::Cross => None,
        Marker::Square => {
            let s = size;
            Some(vec![(-s, -s), (s, -s), (s, s), (-s, s)])
        }
        Marker::Triangle => {
            let h = size * 1.1547;
            Some(vec![(0.0, -h), (-size, h * 0.5), (size, h * 0.5)])
        }
        Marker::TriangleDown => {
            let h = size * 1.1547;
            Some(vec![(0.0, h), (-size, -h * 0.5), (size, -h * 0.5)])
        }
        Marker::Diamond => {
            let s = size * 1.2;
            Some(vec![(0.0, -s), (s, 0.0), (0.0, s), (-s, 0.0)])
        }
        Marker::Plus => {
            let s = size;
            let w = size * 0.3;
            Some(vec![
                (-w, -s), (w, -s), (w, -w), (s, -w), (s, w), (w, w),
                (w, s), (-w, s), (-w, w), (-s, w), (-s, -w), (-w, -w),
            ])
        }
        Marker::Star => {
            let outer = size;
            let inner = size * 0.4;
            let mut pts = Vec::with_capacity(10);
            for i in 0..10 {
                let r = if i % 2 == 0 { outer } else { inner };
                let angle = std::f64::consts::PI * (i as f64) / 5.0 - std::f64::consts::FRAC_PI_2;
                pts.push((r * angle.cos(), r * angle.sin()));
            }
            Some(pts)
        }
        Marker::Pentagon => {
            regular_polygon_points(5, size)
        }
        Marker::Hexagon => {
            regular_polygon_points(6, size)
        }
        Marker::Custom(_) => None, // Custom SVG paths can't be easily rasterized
    }
}

fn regular_polygon_points(sides: usize, radius: f64) -> Option<Vec<(f64, f64)>> {
    let mut pts = Vec::with_capacity(sides);
    let offset = if sides == 5 { -std::f64::consts::FRAC_PI_2 } else { 0.0 };
    for i in 0..sides {
        let angle = std::f64::consts::PI * 2.0 * (i as f64) / (sides as f64) + offset;
        pts.push((radius * angle.cos(), radius * angle.sin()));
    }
    Some(pts)
}

// ─── Backend implementation ──────────────────────────────────────────────────

impl Backend for PngBackend {
    fn draw_line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, style: &LineStyle) {
        let (r, g, b, a) = style.color.to_rgba();
        let a = (a as f64 * style.opacity).round() as u8;

        if let Some(pattern) = style.dash.to_dash_array() {
            self.draw_dashed_line(x1, y1, x2, y2, style.width, &pattern, r, g, b, a);
        } else {
            self.draw_thick_line(x1, y1, x2, y2, style.width, r, g, b, a);
        }
    }

    fn draw_polyline(&mut self, points: &[(f64, f64)], style: &LineStyle) {
        if points.len() < 2 {
            return;
        }
        let (r, g, b, a) = style.color.to_rgba();
        let a = (a as f64 * style.opacity).round() as u8;
        let dash = style.dash.to_dash_array();

        for i in 0..points.len() - 1 {
            let (x0, y0) = points[i];
            let (x1, y1) = points[i + 1];
            if let Some(ref pattern) = dash {
                self.draw_dashed_line(x0, y0, x1, y1, style.width, pattern, r, g, b, a);
            } else {
                self.draw_thick_line(x0, y0, x1, y1, style.width, r, g, b, a);
            }
        }
    }

    fn draw_rect(&mut self, x: f64, y: f64, width: f64, height: f64, style: &FillStyle) {
        let (fr, fg, fb, fa) = style.color.to_rgba();
        let fa = (fa as f64 * style.opacity).round() as u8;
        self.fill_rect_px(x, y, width, height, fr, fg, fb, fa);

        if let Some(ref stroke) = style.stroke {
            let (sr, sg, sb, sa) = stroke.to_rgba();
            self.stroke_rect_px(x, y, width, height, sr, sg, sb, sa, style.stroke_width);
        }
    }

    fn draw_circle(&mut self, cx: f64, cy: f64, r: f64, style: &FillStyle) {
        let (fr, fg, fb, fa) = style.color.to_rgba();
        let fa = (fa as f64 * style.opacity).round() as u8;
        self.fill_circle_px(cx, cy, r, fr, fg, fb, fa);

        if let Some(ref stroke) = style.stroke {
            let (sr, sg, sb, sa) = stroke.to_rgba();
            self.stroke_circle_px(cx, cy, r, sr, sg, sb, sa, style.stroke_width);
        }
    }

    fn draw_text(&mut self, x: f64, y: f64, text: &str, style: &TextStyle) {
        let (r, g, b, a) = style.color.to_rgba();
        let a = (a as f64 * style.opacity).round() as u8;
        let scale = (style.font_size / (FONT_H as f64 + 2.0)).round().max(1.0) as usize;
        self.render_text(x, y, text, scale, &style.anchor, &style.baseline,
                         style.rotation, r, g, b, a);
    }

    fn draw_marker(&mut self, x: f64, y: f64, style: &MarkerStyle) {
        self.render_marker_at(x, y, style);
    }

    fn start_clip(&mut self, x: f64, y: f64, width: f64, height: f64) {
        self.clip_stack.push(ClipRect { x, y, w: width, h: height });
    }

    fn end_clip(&mut self) {
        self.clip_stack.pop();
    }
}

// ─── PNG encoder ─────────────────────────────────────────────────────────────

/// CRC32 lookup table (standard polynomial 0xEDB88320).
fn crc32_table() -> [u32; 256] {
    let mut table = [0u32; 256];
    for i in 0..256 {
        let mut c = i as u32;
        for _ in 0..8 {
            if c & 1 != 0 {
                c = 0xEDB88320 ^ (c >> 1);
            } else {
                c >>= 1;
            }
        }
        table[i] = c;
    }
    table
}

fn crc32(data: &[u8]) -> u32 {
    let table = crc32_table();
    let mut crc = 0xFFFFFFFFu32;
    for &byte in data {
        crc = table[((crc ^ byte as u32) & 0xFF) as usize] ^ (crc >> 8);
    }
    crc ^ 0xFFFFFFFF
}

fn adler32(data: &[u8]) -> u32 {
    let mut a: u32 = 1;
    let mut b: u32 = 0;
    for &byte in data {
        a = (a + byte as u32) % 65521;
        b = (b + a) % 65521;
    }
    (b << 16) | a
}

fn write_chunk(out: &mut Vec<u8>, chunk_type: &[u8; 4], data: &[u8]) {
    let len = data.len() as u32;
    out.extend_from_slice(&len.to_be_bytes());
    out.extend_from_slice(chunk_type);
    out.extend_from_slice(data);
    // CRC over type + data
    let mut crc_data = Vec::with_capacity(4 + data.len());
    crc_data.extend_from_slice(chunk_type);
    crc_data.extend_from_slice(data);
    let crc = crc32(&crc_data);
    out.extend_from_slice(&crc.to_be_bytes());
}

/// Encode RGBA pixel data as a valid PNG file using uncompressed DEFLATE.
fn encode_png_data(width: u32, height: u32, pixels: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();

    // PNG signature
    out.extend_from_slice(&[137, 80, 78, 71, 13, 10, 26, 10]);

    // IHDR chunk
    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&width.to_be_bytes());
    ihdr.extend_from_slice(&height.to_be_bytes());
    ihdr.push(8);  // bit depth
    ihdr.push(6);  // color type: RGBA
    ihdr.push(0);  // compression method
    ihdr.push(0);  // filter method
    ihdr.push(0);  // interlace method
    write_chunk(&mut out, b"IHDR", &ihdr);

    // Prepare raw image data with filter bytes
    let row_bytes = (width * 4) as usize;
    let mut raw = Vec::with_capacity(height as usize * (1 + row_bytes));
    for y in 0..height as usize {
        raw.push(0); // filter type: None
        let start = y * row_bytes;
        raw.extend_from_slice(&pixels[start..start + row_bytes]);
    }

    // Wrap in zlib (uncompressed DEFLATE stored blocks)
    let mut zlib = Vec::new();
    // zlib header: CM=8 (deflate), CINFO=7 (32K window), FCHECK adjusted
    zlib.push(0x78);
    zlib.push(0x01);

    // DEFLATE stored blocks
    let max_block = 65535usize;
    let mut offset = 0;
    while offset < raw.len() {
        let remaining = raw.len() - offset;
        let block_size = remaining.min(max_block);
        let is_final = offset + block_size >= raw.len();
        zlib.push(if is_final { 0x01 } else { 0x00 }); // BFINAL + BTYPE=00
        let len = block_size as u16;
        let nlen = !len;
        zlib.extend_from_slice(&len.to_le_bytes());
        zlib.extend_from_slice(&nlen.to_le_bytes());
        zlib.extend_from_slice(&raw[offset..offset + block_size]);
        offset += block_size;
    }

    // Adler32 of uncompressed data
    let checksum = adler32(&raw);
    zlib.extend_from_slice(&checksum.to_be_bytes());

    // IDAT chunk
    write_chunk(&mut out, b"IDAT", &zlib);

    // IEND chunk
    write_chunk(&mut out, b"IEND", &[]);

    out
}
