// The core of this example was developed by [@bipentihexium](https://github.com/bipentihexium)

use spherical_geometry::{EdgeDirection, Polygon, SphericalPoint};
use std::f32::consts::PI;

const DEFAULT_Y_RANGE: i32 = 50;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let pol = Polygon::new(
        vec![
            SphericalPoint::new(0.0, 0.0),
            SphericalPoint::new(0.0, 0.5),
            SphericalPoint::new(1.2, 0.5),
            SphericalPoint::new(0.8, 0.25),
            SphericalPoint::new(1.2, 0.0),
        ],
        EdgeDirection::CounterClockwise,
    )
    .expect("The polygon should be constructable");
    let y_range = if args.len() < 2 { DEFAULT_Y_RANGE } else { args[1].parse::<i32>().unwrap_or(DEFAULT_Y_RANGE) };
    let x_range = if args.len() < 3 { y_range * 2 } else { args[2].parse::<i32>().unwrap_or(y_range * 2) };
    let (y_start, y_end) = if args.len() < 5 {
        (-PI / 2.0, PI / 2.0)
    } else {
        (args[3].parse::<f32>().unwrap_or(-PI / 2.0), args[4].parse::<f32>().unwrap_or(PI / 2.0))
    };
    let (x_start, x_end) = if args.len() < 7 {
        (-PI, PI)
    } else {
        (args[5].parse::<f32>().unwrap_or(-PI), args[6].parse::<f32>().unwrap_or(PI))
    };
    let ysr = 4;
    let xsr = 2;
    for y in 0..y_range {
        for x in 0..x_range {
            let xfns = (x as f32) / (x_range as f32);
            let yfns = (y as f32) / (y_range as f32);
            let pns = SphericalPoint::new((x_end - x_start) * xfns + x_start, -((y_end - y_start) * yfns + y_start));
            let mut inc = 0u32;
            let samples = xsr * ysr;
            for sy in 0..ysr {
                for sx in 0..xsr {
                    let xf = ((x as f32) + (sx as f32) / (xsr as f32)) / (x_range as f32);
                    let yf = ((y as f32) + (sy as f32) / (ysr as f32)) / (y_range as f32);
                    let p = SphericalPoint::new((x_end - x_start) * xf + x_start, -((y_end - y_start) * yf + y_start));
                    if pol.contains_point(&p).unwrap() {
                        inc += 1;
                    }
                }
            }
            let col = (inc as f32) / (samples as f32);
            let col_r = 255.0 * (2.83 * col * col - 2.36 * col + 0.5);
            let col_g = 255.0 * (2.0 * col - col * col);
            let col_b = 255.0 * 0.7 * (1.0 - col * col);
            let isp = pol.vertices().iter().map(|p2| pns.distance(p2)).fold(f32::INFINITY, |a, b| a.min(b)) < 0.1;
            print!("\x1b[48;2;{};{};{}m{}", col_r as u32, col_g as u32, col_b as u32, if isp { '#' } else { ' ' });
        }
        println!("\x1b[0m");
    }
}
