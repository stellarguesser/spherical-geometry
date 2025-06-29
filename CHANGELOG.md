# Changelog
## 0.4.1 - 2025-06-29 - Internal work and fixes
### ⭐ Added
 - Exposed internal functions for checking the equality of points (mostly a by-product of what can be found in the "Fixed" section)

### 🐛 Fixed
 - Checks for "equal points" are used consistently throughout the crate, which fixes the mismatch between creating a polygon from a list of points and from a text file with serde (the creation could be fine in the former case and break in the latter case before)

## 0.4.0 - 2025-05-25 - Serde support
### ⭐ Added
 - Serde support for structs, optional under a feature flag (https://github.com/stellarguesser/spherical-geometry/pull/11)

## 0.3.0 - 2024-12-19 - More intersections support
### ⭐ Added
 - Checking intersections between great circle arcs (https://github.com/stellarguesser/spherical-geometry/pull/6)
 - Checking intersections between polygons and great circles or great circle arcs (https://github.com/stellarguesser/spherical-geometry/pull/8)
 - Checking if great circles/great circles arcs intersect one another (`true`/`false`, not returning the list of intersections) (https://github.com/stellarguesser/spherical-geometry/pull/8)

### 🐛 Fixed
 - Fixed the example in README (GeCAA Theory task 7) (https://github.com/stellarguesser/spherical-geometry/pull/5)

### 🔧 Improved
- The documentation now uses links to mentioned structs or functions instead of code block references (https://github.com/stellarguesser/spherical-geometry/pull/9)

## 0.2.0 - 2024-10-07 - Polygons support
This release focused on bringing in polygons support, but that required adding several other features :D

### ⭐ Added
 - Polygons construction from vertices
 - Checking if a point is inside a polygon
 - A function to get the closest point on an arc to a given point
 - A function to get the angular distance between two points
 - A function to construct a great circle perpendicular to another great circle or a great circle arc
 - An example showcasing the use of the `Polygon` API
 - The README file now includes an example of using the `Polygon::contains_point` function for determining which stars are inside a constellation.

### 🐛 Fixed
 - Identical great circles are now checked by using the circles' precomputed normals. Before they were checked using new normals, which were however not normalized, leading to wrong results when circles were defined by points close to each other.

### 🔧 Improved
 - The wording of the documentation was changed in several places.

## 0.1.0 - 2024-10-02 - Initial release
This is the initial release of the crate after splitting it away from another codebase.

### ⭐ Added
 - Support for points on the sphere, including spherical ↔ cartesian conversion, (approximate) equality check, and distance between points (metric) functions
 - Support for great circles, including construction from two points and construction from an arc
 - Support for great circle arcs, including construction from two points, checking if it contains a point, getting an intersection with a great circle, and getting a clamped intersection with great circle (returning the closest endpoint if no intersection is on the arc)
