use crate::{GreatCircleArc, SphericalError, SphericalPoint, VEC_LEN_IS_ZERO};

/// Specifies the direction in which an edge is defined.
///
/// Flipping the direction will cause a polygon to be the complement of what you would expect to be your polygon (if you wanted to create a small triangle around the North Pole but flipped the edge orientation, you would define the polygon to be everywhere apart from the North Pole).
///
/// # Determining the direction
/// ## Algorithmic method
/// Imagine you are **inside** the sphere and are looking at the polygon. Now imagine a reference point inside the polygon and find the closest edge to it.
/// Now the direction is the same as the direction of the edge with reference to the point chosen.
///
/// ## More intuitive method
/// Imagine you are standing on the **inside** surface of the sphere, your head pointing in the direction of the centre of the sphere.
/// If you were to walk along the edge and the inside of the polygon was on your left choose `CounterClockwise`, else choose `Clockwise`.
#[derive(Clone, Copy)]
pub enum EdgeDirection {
    Clockwise,
    CounterClockwise,
}

/// A polygon on a unit sphere, given by its vertices and the edge direction
pub struct Polygon {
    vertices: Vec<SphericalPoint>,
    edges_direction: EdgeDirection,
}

impl Polygon {
    /// Creates a new polygon with the vertices and edges direction provided.
    ///
    /// # Important
    /// If the final vertex is not equal to the first one, it will be added automatically.
    ///
    /// Flipping the direction of the edges will cause a polygon to be the complement of what you would expect to be your polygon (if you wanted to create a small triangle around the North Pole but flipped the edge orientation, you would define the polygon to be everywhere apart from the North Pole).
    ///
    /// # Errors
    /// If any edge is defined by essentially equal or antipodal points, returns `SphericalError::AntipodalOrTooClosePoints` as in the case of identical or antipodal points the great circle (and therefore also the edge) is not uniquely defined.
    pub fn new(vertices_in: Vec<SphericalPoint>, edges_direction: EdgeDirection) -> Result<Self, SphericalError> {
        let mut vertices = vertices_in;
        if !vertices[0].approximately_equals(&vertices[vertices.len() - 1], VEC_LEN_IS_ZERO) {
            // The last vertex is not the same as the first one -> insert the first one to the back
            vertices.push(vertices[0]);
        }
        for i in 0..vertices.len() - 1 {
            if vertices[i].cartesian().cross(&vertices[i + 1].cartesian()).magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
                return Err(SphericalError::AntipodalOrTooClosePoints);
            }
        }
        Ok(Self { vertices, edges_direction })
    }

    /// Returns a reference to the vertices list
    pub fn vertices(&self) -> &Vec<SphericalPoint> {
        &self.vertices
    }

    /// Returns the edges direction
    pub fn edges_direction(&self) -> EdgeDirection {
        self.edges_direction
    }

    ///
    ///
    /// # Errors
    /// If any of the edges fails to be constructed as a `GreatCircleArc`, returns the corresponding error. This should however never happen, as that is checked when the polygon is constructed.
    pub fn contains_point(&self, point: &SphericalPoint) -> Result<bool, SphericalError> {
        // Algorithm description:
        // 1) Find the closest edge by finding an intersection with each of the edges with a great circle perpendicular to it. Use the clamped intersection, returning one of the endpoints in case of a miss.
        // 2) Determine if the closest edge is in the correct orientation

        // Step 1
        let mut closest_edge_i = 0;
        let mut closest_edge_dist_metric = f32::INFINITY;
        for i in 0..self.vertices.len() - 1 {
            let edge = GreatCircleArc::new(self.vertices[i], self.vertices[i + 1])?;
            if edge.contains_point(point) {
                return Ok(true);
            }
            let edge_distance_metric = match edge.perpendicular_circle_through_point(point) {
                Ok(circle) => {
                    let closest_intersection = edge
                        .intersect_great_circle_clamped(&circle)
                        .expect("Perpendicular great circle must not be identical to the original arc.");
                    if closest_intersection.is_empty() {
                        continue;
                    }
                    closest_intersection[0].minus_cotan_distance(point)
                }
                Err(SphericalError::AntipodalOrTooClosePoints) => {
                    // The point is essentially the pole of the arc, so it is basically PI/2 radians away -> distance metric = -1/tan(PI/2) = 0
                    0.0
                }
                Err(err) => return Err(err),
            };
            if edge_distance_metric < closest_edge_dist_metric {
                closest_edge_i = i;
                closest_edge_dist_metric = edge_distance_metric;
            }
        }

        // Step 2
        let closest_edge_normal = self.vertices[closest_edge_i].cartesian().cross(&self.vertices[closest_edge_i + 1].cartesian());
        let cos_angle = closest_edge_normal.dot(&point.cartesian());
        let is_inside = match self.edges_direction {
            EdgeDirection::Clockwise => cos_angle >= 0.0,
            EdgeDirection::CounterClockwise => cos_angle <= 0.0,
        };

        Ok(is_inside)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn is_point_inside() {
        let polygon_1 = Polygon::new(
            vec![
                SphericalPoint::new(0.0, PI / 3.0),
                SphericalPoint::new(-2.0 * PI / 3.0, PI / 3.0),
                SphericalPoint::new(2.0 * PI / 3.0, PI / 3.0),
            ],
            EdgeDirection::CounterClockwise,
        )
        .expect("The polygon should be constructable");
        let north_pole = SphericalPoint::new(0.0, PI / 2.0);
        assert!(polygon_1.contains_point(&north_pole).expect("It should be possible to determine if the point is inside the polygon"));
        let south_pole = SphericalPoint::new(0.0, -PI / 2.0);
        assert!(!polygon_1.contains_point(&south_pole).expect("It should be possible to determine if the point is inside the polygon"));
    }
}
