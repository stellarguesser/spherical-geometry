# Changelog
## 0.1.0 - 2024-10-02 - Initial release
This is the initial release of the crate after splitting it away from another codebase.

### ⭐ Added
 - Support for points on the sphere, including spherical ↔ cartesian conversion, (approximate) equality check, and distance between points (metric) functions
 - Support for great circles, including construction from two points and construction from an arc
 - Support for great circle arcs, including construction from two points, checking if it contains a point, getting an intersection with a great circle, and getting a clamped intersection with great circle (returning the closest endpoint if no intersection is on the arc)