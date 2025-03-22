//! # Spherical Geometry
//!
//! `spherical_geometry` is a library for handling geometry on the surface of a sphere.
//!
//! This library combines spherical and vector geometry to perform operations with points, [great circles](https://en.wikipedia.org/wiki/Great_circle), great circle arcs... A great circle is an equivalent of a straight line in planar geometry - it is the shortest path between two points on a sphere.
//!
//! Read more about it at [https://github.com/stellarguesser/spherical-geometry](https://github.com/stellarguesser/spherical-geometry)
//!
//! If it is not obvious how to use a function and there is not an example in the documentation, please open an issue/pull request in the repository. In the meantime before it gets added, you can try to take a look at the unit tests, which usually show how to use a function, which can serve as an example.
//!
//! ## Feature flags
//! - `serde`: Enables serialization and deserialization of types using Serde.

pub mod great_circle;
pub mod great_circle_arc;
pub mod point;
pub mod polygon;

pub use great_circle::GreatCircle;
pub use great_circle_arc::GreatCircleArc;
pub use point::SphericalPoint;
pub use polygon::{EdgeDirection, Polygon};

pub(crate) const IDENTICAL_POINTS: f32 = 10e-5;
pub(crate) const VEC_LEN_IS_ZERO: f32 = 10e-6;

#[derive(Debug)]
pub enum SphericalError {
    AntipodalOrTooClosePoints,
    IdenticalGreatCircles,
    PoleAndPointNotNormal,
}
