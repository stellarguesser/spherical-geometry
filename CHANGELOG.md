# Changelog
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

### Improved
 - The wording of the documentation was changed in several places.

## 0.1.0 - 2024-10-02 - Initial release
This is the initial release of the crate after splitting it away from another codebase.

### ⭐ Added
 - Support for points on the sphere, including spherical ↔ cartesian conversion, (approximate) equality check, and distance between points (metric) functions
 - Support for great circles, including construction from two points and construction from an arc
 - Support for great circle arcs, including construction from two points, checking if it contains a point, getting an intersection with a great circle, and getting a clamped intersection with great circle (returning the closest endpoint if no intersection is on the arc)
