# Examples
This folder contains examples for using the library. You can run each of them with `cargo run --example <example-name>` (for example `cargo run --example cli_polygon_test`).

## cli polygon test
The core of this example was developed by [@bipentihexium](https://github.com/bipentihexium). It is a cli utility to test points on the sphere if they are inside a polygon. The default polygon is a non-convex one that kept causing issues in development :D It accepts several command line arguments, in order:
 - `<height>` - how many rows to print the result into (kind of resolution). If it can not be parsed to an `i32`, defaults to `50`.
 - `<width>` - how many columns to print the result into (kind of resolution). If it can not be parsed to an `i32`, defaults to `height*2`.
 - `<dec-start>` - the declination to start at, in radians. If it can not be parsed to an `f32`, defaults to `-PI/2.0`.
 - `<dec-end>` - the declination to end at, in radians. If it can not be parsed to an `f32`, defaults to `PI/2.0`.
 - `<ra-start>` - the right ascension to start at, in radians. If it can not be parsed to an `f32`, defaults to `-PI`.
 - `<ra-end>` - the right ascension to end at, in radians. If it can not be parsed to an `f32`, defaults to `PI`.

Points considered to be inside the polygon are yellow, those outside are purple, and colours in between indicate a state in between - it uses MSAA, so each printed point actually checks multiple points in its area and averages the result.