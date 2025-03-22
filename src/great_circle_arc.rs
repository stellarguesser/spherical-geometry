use crate::{GreatCircle, SphericalError, SphericalPoint, VEC_LEN_IS_ZERO};
use nalgebra::Vector3;
#[cfg(feature = "serde")]
use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// A great circle on a unit sphere, given by two points on it
#[derive(Clone, Copy, Debug)]
pub struct GreatCircleArc {
    start: SphericalPoint,
    end: SphericalPoint,
    normal: Vector3<f32>,
}

impl GreatCircleArc {
    /// Creates a new great circle arc passing through the two points provided, taking the shorter of the two possible paths
    ///
    /// # Errors
    /// If the points are essentially equal or essentially antipodal, returns [SphericalError::AntipodalOrTooClosePoints] as in the case of identical or antipodal points the great circle (and therefore also the arc) is not uniquely defined
    pub fn new(point1: SphericalPoint, point2: SphericalPoint) -> Result<Self, SphericalError> {
        if point1.cartesian().cross(&point2.cartesian()).magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
            return Err(SphericalError::AntipodalOrTooClosePoints);
        }
        Ok(Self {
            start: point1,
            end: point2,
            normal: point1.cartesian().cross(&point2.cartesian()).normalize(),
        })
    }

    /// Returns the start of the great circle arc
    pub fn start(&self) -> SphericalPoint {
        self.start
    }

    /// Returns the end of the great circle arc
    pub fn end(&self) -> SphericalPoint {
        self.end
    }

    /// Returns the vector normal to the great circle containing the arc
    pub fn normal(&self) -> Vector3<f32> {
        self.normal
    }

    /// Checks if the great circle arc contains the provided point
    pub fn contains_point(&self, point: &SphericalPoint) -> bool {
        let tolerance: f32 = 10e-5;
        let great_circle = GreatCircle::from_arc(self);
        if !great_circle.contains_point(point) {
            return false;
        }
        // If the point is approximately equal to either of the ends, it obviously is on the arc
        if self.start.approximately_equals(point, tolerance) || self.end.approximately_equals(point, tolerance) {
            return true;
        }
        // If angle AOP + angle POB = angle AOB, then the point is on the arc -> Check for cos(AOP + POB) = cos(AOP) * cos(POB) - sin(AOP) * sin(POB)
        // This way we avoid relatively costly inverse trigonometric functions
        let cos_aob = self.start.cartesian().dot(&self.end.cartesian());
        let cos_aop = self.start.cartesian().dot(&point.cartesian());
        let cos_pob = point.cartesian().dot(&self.end.cartesian());

        #[cfg(test)]
        dbg!(cos_aob, cos_aop, cos_pob);

        // If either of the cosines is smaller than cos(angle AOB), then the angular distance from one of the endings of the arc to the point is greater than the length of the arc -> the point must be outside
        // This avoids the issue of points that are opposite the arc
        if cos_aob - cos_aop > tolerance || cos_aob - cos_pob > tolerance {
            return false;
        }
        let sin_aop = self.start.cartesian().cross(&point.cartesian()).magnitude();
        let sin_pob = point.cartesian().cross(&self.end.cartesian()).magnitude();
        let cos_aob_calc = cos_aop * cos_pob - sin_aop * sin_pob;

        #[cfg(test)]
        dbg!(cos_aob_calc);

        (cos_aob - cos_aob_calc).abs() < tolerance
    }

    /// Creates a new great circle passing through the provided point and perpendicular to the current circle arc
    ///
    /// # Errors
    /// If the point and the pole of the current circle arc are essentially equal or essentially antipodal, returns [SphericalError::AntipodalOrTooClosePoints] as in the case of identical or antipodal points the great circle is not uniquely defined
    pub fn perpendicular_circle_through_point(&self, point: &SphericalPoint) -> Result<GreatCircle, SphericalError> {
        let point_1 = SphericalPoint::from_cartesian_vector3(self.normal());
        GreatCircle::new(point_1, *point)
    }

    /// Returns the intersections of this great circle arc with a great circle
    ///
    /// # Errors
    /// If the great circle and the great circle containing the arc are (essentially) parallel (equal to each other), returns [SphericalError::IdenticalGreatCircles] as then there is an infinite amount of intersections. You can handle this error as an equivalent of "all points on the arc are intersections".
    pub fn intersect_great_circle(&self, other: &GreatCircle) -> Result<Vec<SphericalPoint>, SphericalError> {
        let normal1 = self.normal();
        let normal2 = other.normal();

        let res = normal1.cross(&normal2);
        if res.magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
            return Err(SphericalError::IdenticalGreatCircles);
        }
        let res_norm = res.normalize();
        let point_1 = SphericalPoint::from_cartesian_vector3(res_norm);
        let point_2 = SphericalPoint::from_cartesian_vector3(-res_norm);

        #[cfg(test)]
        dbg!(point_1);
        #[cfg(test)]
        dbg!(point_2);

        let mut intersections = Vec::new();
        if self.contains_point(&point_1) {
            intersections.push(point_1);
        }
        if self.contains_point(&point_2) {
            intersections.push(point_2);
        }
        Ok(intersections)
    }

    /// Checks if there exists an intersection between this great circle arc and the provided great circle
    ///
    /// # Errors
    /// Only propagates errors originating from [Self::intersect_great_circle], but handles [SphericalError::IdenticalGreatCircles] internally
    pub fn intersects_great_circle(&self, circle: &GreatCircle) -> Result<bool, SphericalError> {
        match self.intersect_great_circle(circle) {
            Ok(intersections) => Ok(!intersections.is_empty()),
            Err(err) => {
                match err {
                    SphericalError::IdenticalGreatCircles => {
                        // They are parallel -> the arc is a part of the circle
                        Ok(true)
                    }
                    _ => Err(err),
                }
            }
        }
    }

    /// Returns the intersections of this great circle arc with another one
    ///
    /// # Errors
    /// If the great circles containing the arcs are (essentially) parallel (equal to each other) and overlapping, returns [SphericalError::IdenticalGreatCircles] as then there is an infinite amount of intersections.
    ///
    /// Propagates the rest of errors originating from [Self::intersect_great_circle]
    pub fn intersect_great_circle_arc(&self, other: &Self) -> Result<Vec<SphericalPoint>, SphericalError> {
        let circle = GreatCircle::from_arc(other);
        match self.intersect_great_circle(&circle) {
            Ok(intersections) => Ok(intersections.into_iter().filter(|intersection| other.contains_point(intersection)).collect()),
            Err(err) => {
                match err {
                    SphericalError::IdenticalGreatCircles => {
                        if self.contains_point(&other.start) || self.contains_point(&other.end) || other.contains_point(&self.start) || other.contains_point(&self.end) {
                            // They are parallel and overlap each other
                            Err(SphericalError::IdenticalGreatCircles)
                        } else {
                            // They are parallel, but do not overlap each other
                            Ok(Vec::new())
                        }
                    }
                    _ => Err(err),
                }
            }
        }
    }

    /// Checks if there exists an intersection of this arc with the provided one
    ///
    /// # Errors
    /// Propagates errors originating from [Self::intersect_great_circle] apart from [SphericalError::IdenticalGreatCircles] for which the meaning in this context can be determined
    pub fn intersects_great_circle_arc(&self, other: &Self) -> Result<bool, SphericalError> {
        let circle = GreatCircle::from_arc(other);
        match self.intersect_great_circle(&circle) {
            Ok(intersections) => Ok(intersections.into_iter().any(|intersection| other.contains_point(&intersection))),
            Err(err) => {
                match err {
                    SphericalError::IdenticalGreatCircles => {
                        if self.contains_point(&other.start) || self.contains_point(&other.end) || other.contains_point(&self.start) || other.contains_point(&self.end) {
                            // They are parallel and overlap each other
                            Ok(true)
                        } else {
                            // They are parallel, but do not overlap each other
                            Ok(false)
                        }
                    }
                    _ => Err(err),
                }
            }
        }
    }

    /// Returns a point that is closest to the given point using a provided perpendicular circle. This is to be used when one has already constructed the perpendicular circle.
    ///
    /// The provided circle is intended to be perpendicular to the arc as that is the only time this function will return meaningful results. It does not, however, rely on this being the case, so you can use it with any circle if you find it useful.
    ///
    /// # Errors
    /// If the perpendicular great circle and the great circle containing the arc are (essentially) parallel (equal to each other), returns [SphericalError::IdenticalGreatCircles] as then there is an infinite amount of intersections. This should, however, never happen and would indicate a bug in the implementation.
    pub fn closest_point_to_point_with_circle(&self, perpendicular_circle: &GreatCircle, point: &SphericalPoint) -> Result<SphericalPoint, SphericalError> {
        let normal1 = self.normal();
        let normal2 = perpendicular_circle.normal();

        let res = normal1.cross(&normal2);
        if res.magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
            return Err(SphericalError::IdenticalGreatCircles);
        }
        let res_norm = res.normalize();
        let point_1 = SphericalPoint::from_cartesian_vector3(res_norm);
        let point_2 = SphericalPoint::from_cartesian_vector3(-res_norm);

        #[cfg(test)]
        dbg!(point_1);
        #[cfg(test)]
        dbg!(point_2);

        let mut points_considered = vec![self.start, self.end];
        if self.contains_point(&point_1) {
            points_considered.push(point_1);
        }
        if self.contains_point(&point_2) {
            points_considered.push(point_2);
        }

        let mut closest_metric = f32::INFINITY;
        let mut closest_point = points_considered[0];

        for point_c in points_considered {
            let distance_metric = point_c.minus_cotan_distance(point);
            if distance_metric < closest_metric {
                closest_metric = distance_metric;
                closest_point = point_c;
            }
        }
        Ok(closest_point)
    }

    /// Returns the closest point to the arc from a given point.
    ///
    /// # Errors
    /// If the perpendicular great circle can not be constructed (usually because the given point is a pole of the circle containing the arc), returns [SphericalError::AntipodalOrTooClosePoints] as in the case of identical or antipodal points the great circle is not uniquely defined.
    ///
    /// If the perpendicular great circle and the great circle containing the arc are (essentially) parallel (equal to each other), returns [SphericalError::IdenticalGreatCircles] as then there is an infinite amount of intersections. This should, however, never happen and would indicate a bug in the implementation.
    pub fn closest_point_to_point(&self, point: &SphericalPoint) -> Result<SphericalPoint, SphericalError> {
        let perpendicular_circle = self.perpendicular_circle_through_point(point)?;

        let normal1 = self.normal();
        let normal2 = perpendicular_circle.normal();

        let res = normal1.cross(&normal2);
        if res.magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
            return Err(SphericalError::IdenticalGreatCircles);
        }
        let res_norm = res.normalize();
        let point_1 = SphericalPoint::from_cartesian_vector3(res_norm);
        let point_2 = SphericalPoint::from_cartesian_vector3(-res_norm);

        #[cfg(test)]
        dbg!(point_1);
        #[cfg(test)]
        dbg!(point_2);

        let mut points_considered = vec![self.start, self.end];
        if self.contains_point(&point_1) {
            points_considered.push(point_1);
        }
        if self.contains_point(&point_2) {
            points_considered.push(point_2);
        }

        let mut closest_metric = f32::INFINITY;
        let mut closest_point = points_considered[0];

        for point_c in points_considered {
            let distance_metric = point_c.minus_cotan_distance(point);
            if distance_metric < closest_metric {
                closest_metric = distance_metric;
                closest_point = point_c;
            }
        }
        Ok(closest_point)
    }

    /// Returns the intersections of the arc with the great circle, clamped to the arc. If there are no intersections, the endpoint closest to the potential intersections (of the great circle and the arc extended into a great circle) is returned.
    ///
    /// # Errors
    /// If the great circle and the great circle containing the arc are (essentially) parallel (equal to each other), returns [SphericalError::IdenticalGreatCircles] as then there is an infinite amount of intersections. You can handle this error as an equivalent of "all points on the arc are intersections".
    pub fn intersect_great_circle_clamped(&self, circle: &GreatCircle) -> Result<Vec<SphericalPoint>, SphericalError> {
        let normal1 = self.normal();
        let normal2 = circle.normal();

        let res = normal1.cross(&normal2);
        if res.magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
            return Err(SphericalError::IdenticalGreatCircles);
        }
        let res_norm = res.normalize();
        let point_1 = SphericalPoint::from_cartesian_vector3(res_norm);
        let point_2 = SphericalPoint::from_cartesian_vector3(-res_norm);

        #[cfg(test)]
        dbg!(point_1);
        #[cfg(test)]
        dbg!(point_2);

        let mut intersections = Vec::new();
        if self.contains_point(&point_1) {
            intersections.push(point_1);
        }
        if self.contains_point(&point_2) {
            intersections.push(point_2);
        }
        if intersections.is_empty() {
            let start_1_distance = self.start.minus_cotan_distance(&point_1);
            let start_2_distance = self.start.minus_cotan_distance(&point_2);
            let start_distance = start_1_distance.min(start_2_distance);
            let end_1_distance = self.end.minus_cotan_distance(&point_1);
            if end_1_distance < start_distance {
                intersections.push(self.end);
            } else {
                let end_2_distance = self.end.minus_cotan_distance(&point_2);
                if end_2_distance < start_distance {
                    intersections.push(self.end);
                } else {
                    intersections.push(self.start);
                }
            }
        }
        Ok(intersections)
    }

    /// Returns the intersections of the arc with the great circle, clamped to the arc. If there are no intersections, the endpoint closest to the point provided is returned.
    ///
    /// # Errors
    /// If the great circle and the great circle containing the arc are (essentially) parallel (equal to each other), returns [SphericalError::IdenticalGreatCircles] as then there is an infinite amount of intersections. You can handle this error as an equivalent of "all points on the arc are intersections".
    pub fn intersect_great_circle_clamped_closest_to_point(&self, circle: &GreatCircle, point: &SphericalPoint) -> Result<Vec<SphericalPoint>, SphericalError> {
        let normal1 = self.normal();
        let normal2 = circle.normal();

        let res = normal1.cross(&normal2);
        if res.magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
            return Err(SphericalError::IdenticalGreatCircles);
        }
        let res_norm = res.normalize();
        let point_1 = SphericalPoint::from_cartesian_vector3(res_norm);
        let point_2 = SphericalPoint::from_cartesian_vector3(-res_norm);

        #[cfg(test)]
        dbg!(point_1);
        #[cfg(test)]
        dbg!(point_2);

        let mut intersections = Vec::new();
        if self.contains_point(&point_1) {
            intersections.push(point_1);
        }
        if self.contains_point(&point_2) {
            intersections.push(point_2);
        }
        if intersections.is_empty() {
            let start_distance = self.start.minus_cotan_distance(point);
            let end_distance = self.end.minus_cotan_distance(point);

            if end_distance < start_distance {
                intersections.push(self.end);
            } else {
                intersections.push(self.start);
            }
        }
        Ok(intersections)
    }
}

#[cfg(feature = "serde")]
impl Serialize for GreatCircleArc {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let data = (self.start, self.end);
        data.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for GreatCircleArc {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let (start, end) = <(SphericalPoint, SphericalPoint)>::deserialize(deserializer)?;
        GreatCircleArc::new(start, end).map_err(|e| serde::de::Error::custom(format!("invalid great circle arc: {:?}", e)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn contains_point() {
        let equator_north_to_west = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, 0.0)).expect("The points are far enough");
        let north = SphericalPoint::new(0.0, 0.0);
        let northwest = SphericalPoint::new(PI / 4.0, 0.0);
        let west = SphericalPoint::new(PI / 2.0, 0.0);
        let southeast = SphericalPoint::new(-3.0 * PI / 4.0, 0.0);
        let outside_opposite = SphericalPoint::new(-3.0 * PI / 4.0 - PI / 10.0, 0.0);
        let outside_in_plane = SphericalPoint::new(-PI / 4.0, 0.0);
        let outside_above = SphericalPoint::new(PI / 4.0, PI / 7.0);
        let outside_total = SphericalPoint::new(PI, -PI / 3.0);

        assert!(equator_north_to_west.contains_point(&north));
        assert!(equator_north_to_west.contains_point(&northwest));
        assert!(equator_north_to_west.contains_point(&west));
        assert!(!equator_north_to_west.contains_point(&outside_opposite));
        assert!(!equator_north_to_west.contains_point(&southeast));
        assert!(!equator_north_to_west.contains_point(&outside_in_plane));
        assert!(!equator_north_to_west.contains_point(&outside_above));
        assert!(!equator_north_to_west.contains_point(&outside_total));

        let tolerance = 10e-5;
        for i in 0..360 {
            let angle = 2.0 * PI / 360.0 * (i as f32);
            let point = SphericalPoint::new(angle, 0.0);
            dbg!(point);
            assert_eq!(equator_north_to_west.contains_point(&point), angle < PI / 2.0 || (PI / 2.0 - angle).abs() < tolerance);
        }
    }

    #[test]
    fn intersect_great_circle() {
        let tolerance = 10e-4;
        let arc_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let circle_1 = GreatCircle::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_1 = arc_1.intersect_great_circle(&circle_1).expect("The circles are not parallel");
        assert_eq!(intersections_1.len(), 1);
        assert!(intersections_1[0].approximately_equals(&SphericalPoint::new(0.0, 0.0), tolerance));

        let arc_2 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 6.0), SphericalPoint::new(0.0, PI / 4.0)).expect("The points are far enough");
        let circle_2 = GreatCircle::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_2 = arc_2.intersect_great_circle(&circle_2).expect("The circles are not parallel");
        assert!(intersections_2.is_empty());

        let arc_3 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, PI / 4.0)).expect("The points are far enough");
        let circle_3 = GreatCircle::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_3 = arc_3.intersect_great_circle(&circle_3).expect("The circles are not parallel");
        assert_eq!(intersections_3.len(), 1);
        assert!(intersections_3[0].approximately_equals(&SphericalPoint::new(0.0, 0.0), tolerance));

        let arc_4 = GreatCircleArc::new(SphericalPoint::new(PI / 5.0, -PI / 7.0), SphericalPoint::new(PI / 2.0, PI / 4.0)).expect("The points are far enough");
        let circle_4 = GreatCircle::new(SphericalPoint::new(PI / 5.0, PI / 4.0), SphericalPoint::new(PI / 5.0, -PI / 4.0)).expect("The points are far enough");
        let intersections_4 = arc_4.intersect_great_circle(&circle_4).expect("The circles are not parallel");
        assert_eq!(intersections_4.len(), 1);
        assert!(intersections_4[0].approximately_equals(&SphericalPoint::new(PI / 5.0, -PI / 7.0), tolerance));
    }

    #[test]
    fn intersect_great_circle_arc() {
        let tolerance = 10e-4;

        let arc_1_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let arc_1_2 = GreatCircleArc::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_1 = arc_1_1.intersect_great_circle_arc(&arc_1_2).expect("The circles are not parallel");
        assert_eq!(intersections_1.len(), 1);
        assert!(intersections_1[0].approximately_equals(&SphericalPoint::new(0.0, 0.0), tolerance));

        let arc_2_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 6.0), SphericalPoint::new(0.0, PI / 4.0)).expect("The points are far enough");
        let arc_2_2 = GreatCircleArc::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_2 = arc_2_1.intersect_great_circle_arc(&arc_2_2).expect("The circles are not parallel");
        assert!(intersections_2.is_empty());

        let arc_3_1 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, PI / 4.0)).expect("The points are far enough");
        let arc_3_2 = GreatCircleArc::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_3 = arc_3_1.intersect_great_circle_arc(&arc_3_2).expect("The circles are not parallel");
        assert_eq!(intersections_3.len(), 1);
        assert!(intersections_3[0].approximately_equals(&SphericalPoint::new(0.0, 0.0), tolerance));

        // Arcs miss but the great circles do not
        let arc_4_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let arc_4_2 = GreatCircleArc::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(PI / 3.0, 0.0)).expect("The points are far enough");
        let intersections_4 = arc_4_1.intersect_great_circle_arc(&arc_4_2).expect("The circles are not parallel");
        assert!(intersections_4.is_empty());

        // The arcs are on the opposite sides of the sphere
        let arc_5_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let arc_5_2 = GreatCircleArc::new(SphericalPoint::new(3.0 * PI / 4.0, 0.0), SphericalPoint::new(5.0 * PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_5 = arc_5_1.intersect_great_circle_arc(&arc_5_2).expect("The circles are not parallel");
        assert!(intersections_5.is_empty());

        // Parallel, but not overlapping
        let arc_6_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let arc_6_2 = GreatCircleArc::new(SphericalPoint::new(PI, PI / 4.0), SphericalPoint::new(PI, -PI / 4.0)).expect("The points are far enough");
        let intersections_6 = arc_6_1.intersect_great_circle_arc(&arc_6_2).expect("The arcs do not overlap");
        assert!(intersections_6.is_empty());

        // Parallel and overlapping, but none contains the other
        let arc_7_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let arc_7_2 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 3.0), SphericalPoint::new(0.0, -PI / 5.0)).expect("The points are far enough");
        let intersections_7 = arc_7_1.intersect_great_circle_arc(&arc_7_2);
        assert!(matches!(intersections_7, Err(SphericalError::IdenticalGreatCircles)));

        // Parallel and overlapping, but first contains the second
        let arc_8_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let arc_8_2 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 6.0), SphericalPoint::new(0.0, -PI / 5.0)).expect("The points are far enough");
        let intersections_8 = arc_8_1.intersect_great_circle_arc(&arc_8_2);
        assert!(matches!(intersections_8, Err(SphericalError::IdenticalGreatCircles)));

        // Parallel and overlapping, but second contains the first
        let arc_9_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 5.0)).expect("The points are far enough");
        let arc_9_2 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 3.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let intersections_8 = arc_9_1.intersect_great_circle_arc(&arc_9_2);
        assert!(matches!(intersections_8, Err(SphericalError::IdenticalGreatCircles)));
    }

    #[test]
    fn intersect_great_circle_clamped() {
        let tolerance = 10e-4;
        let arc_1 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 4.0), SphericalPoint::new(0.0, -PI / 4.0)).expect("The points are far enough");
        let circle_1 = GreatCircle::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_1 = arc_1.intersect_great_circle_clamped(&circle_1).expect("The circles are not parallel");
        assert_eq!(intersections_1.len(), 1);
        assert!(intersections_1[0].approximately_equals(&SphericalPoint::new(0.0, 0.0), tolerance));

        let arc_2 = GreatCircleArc::new(SphericalPoint::new(0.0, PI / 6.0), SphericalPoint::new(0.0, PI / 4.0)).expect("The points are far enough");
        let circle_2 = GreatCircle::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_2 = arc_2.intersect_great_circle_clamped(&circle_2).expect("The circles are not parallel");
        assert_eq!(intersections_2.len(), 1);
        assert!(intersections_2[0].approximately_equals(&SphericalPoint::new(0.0, PI / 6.0), tolerance));

        let arc_3 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, PI / 4.0)).expect("The points are far enough");
        let circle_3 = GreatCircle::new(SphericalPoint::new(PI / 4.0, 0.0), SphericalPoint::new(-PI / 4.0, 0.0)).expect("The points are far enough");
        let intersections_3 = arc_3.intersect_great_circle_clamped(&circle_3).expect("The circles are not parallel");
        assert_eq!(intersections_3.len(), 1);
        assert!(intersections_3[0].approximately_equals(&SphericalPoint::new(0.0, 0.0), tolerance));

        let arc_4 = GreatCircleArc::new(SphericalPoint::new(PI / 5.0, -PI / 7.0), SphericalPoint::new(PI / 2.0, PI / 4.0)).expect("The points are far enough");
        let circle_4 = GreatCircle::new(SphericalPoint::new(PI / 5.0, PI / 4.0), SphericalPoint::new(PI / 5.0, -PI / 4.0)).expect("The points are far enough");
        let intersections_4 = arc_4.intersect_great_circle_clamped(&circle_4).expect("The circles are not parallel");
        assert_eq!(intersections_4.len(), 1);
        assert!(intersections_4[0].approximately_equals(&SphericalPoint::new(PI / 5.0, -PI / 7.0), tolerance));

        let arc_5 = GreatCircleArc::new(SphericalPoint::new(PI / 5.0, PI / 7.0), SphericalPoint::new(PI / 2.0, PI / 6.0)).expect("The points are far enough");
        let circle_5 = GreatCircle::new(SphericalPoint::new(PI / 2.0 + PI / 5.0, 0.0), SphericalPoint::new(PI / 2.0 + PI / 5.0, PI / 4.0)).expect("The points are far enough");
        let intersections_5 = arc_5.intersect_great_circle_clamped(&circle_5).expect("The circles are not parallel");
        assert_eq!(intersections_5.len(), 1);
        assert!(intersections_5[0].approximately_equals(&SphericalPoint::new(PI / 2.0, PI / 6.0), tolerance));
    }

    #[test]
    fn perpendicular_circle_through_point() {
        let tolerance = 10e-6;

        let equator = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(-PI / 2.0, 0.0)).expect("The points are fairly far away");
        let greenwich_meridian = equator.perpendicular_circle_through_point(&SphericalPoint::new(0.0, 0.0)).expect("The points are fairly far away");
        let new_pole = SphericalPoint::from_cartesian_vector3(greenwich_meridian.normal());
        let new_pole_corr_1 = SphericalPoint::new(PI / 2.0, 0.0);
        let new_pole_corr_2 = SphericalPoint::new(-PI / 2.0, 0.0);
        assert!(new_pole.approximately_equals(&new_pole_corr_1, tolerance) || new_pole.approximately_equals(&new_pole_corr_2, tolerance));

        let equator = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(-PI / 2.0, 0.0)).expect("The points are fairly far away");
        let greenwich_meridian = equator.perpendicular_circle_through_point(&SphericalPoint::new(0.0, PI / 2.0)); // Telling it to pass through the North Pole, which it should always do
        assert!(greenwich_meridian.is_err());

        let circle_1 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, PI / 6.0)).expect("The points are fairly far away");
        let circle_2 = circle_1
            .perpendicular_circle_through_point(&SphericalPoint::new(0.0, PI / 2.0))
            .expect("The points are fairly far away");
        let new_pole = SphericalPoint::from_cartesian_vector3(circle_2.normal());
        let new_pole_corr_1 = SphericalPoint::new(0.0, 0.0);
        let new_pole_corr_2 = SphericalPoint::new(PI, 0.0);
        assert!(new_pole.approximately_equals(&new_pole_corr_1, tolerance) || new_pole.approximately_equals(&new_pole_corr_2, tolerance));

        let circle_1 = GreatCircleArc::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, PI / 6.0)).expect("The points are fairly far away");
        let circle_2 = circle_1.perpendicular_circle_through_point(&SphericalPoint::new(0.0, 0.0)).expect("The points are fairly far away");
        let new_pole = SphericalPoint::from_cartesian_vector3(circle_2.normal());
        let new_pole_corr_1 = SphericalPoint::new(PI / 2.0, PI / 6.0);
        let new_pole_corr_2 = SphericalPoint::new(-PI / 2.0, -PI / 6.0);
        assert!(new_pole.approximately_equals(&new_pole_corr_1, tolerance) || new_pole.approximately_equals(&new_pole_corr_2, tolerance));
    }

    #[cfg(feature = "serde")]
    mod serde_tests {
        use super::*;
        use serde_json;

        #[test]
        fn test_serde() {
            let point1 = SphericalPoint::new(1.0, 0.5);
            let point2 = SphericalPoint::new(2.0, 1.0);
            let orig = GreatCircleArc::new(point1, point2).expect("Valid points should create a GreatCircleArc");
            let ser = serde_json::to_string(&orig).expect("Serialization failed");
            let deser: GreatCircleArc = serde_json::from_str(&ser).expect("Deserialization failed");

            assert!((orig.start().ra() - deser.start().ra()).abs() < f32::EPSILON, "Start RA values do not match");
            assert!((orig.start().dec() - deser.start().dec()).abs() < f32::EPSILON, "Start Dec values do not match");
            assert!((orig.end().ra() - deser.end().ra()).abs() < f32::EPSILON, "End RA values do not match");
            assert!((orig.end().dec() - deser.end().dec()).abs() < f32::EPSILON, "End Dec values do not match");
            assert!((orig.normal() - deser.normal()).magnitude() < f32::EPSILON, "Normal vector does not match after deserialization");
        }
    }
}
