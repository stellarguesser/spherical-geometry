use crate::{GreatCircleArc, SphericalError, SphericalPoint, VEC_LEN_IS_ZERO};
use nalgebra::Vector3;
#[cfg(feature = "serde")]
use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// A great circle on a unit sphere, given by two points on it
#[derive(Clone, Copy, Debug)]
pub struct GreatCircle {
    start: SphericalPoint,
    end: SphericalPoint,
    normal: Vector3<f32>,
}

impl GreatCircle {
    /// Creates a new great circle passing through the two points provided
    ///
    /// # Errors
    /// If the points are essentially equal or essentially antipodal, returns [SphericalError::AntipodalOrTooClosePoints] as in the case of identical or antipodal points the great circle is not uniquely defined
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

    /// Creates a great circle from an arc
    pub fn from_arc(arc: &GreatCircleArc) -> Self {
        Self {
            start: arc.start(),
            end: arc.end(),
            normal: arc.normal(),
        }
    }

    /// Creates a new great circle passing through the provided point and perpendicular to the current circle
    ///
    /// # Errors
    /// If the point and the pole of the current circle are essentially equal or essentially antipodal, returns [SphericalError::AntipodalOrTooClosePoints] as in the case of identical or antipodal points the great circle is not uniquely defined.
    pub fn perpendicular_through_point(&self, point: &SphericalPoint) -> Result<Self, SphericalError> {
        let point_1 = SphericalPoint::from_cartesian_vector3(self.normal());
        Self::new(point_1, *point)
    }

    /// Returns one of the points defining the great circle
    pub fn start(&self) -> SphericalPoint {
        self.start
    }

    /// Returns the other of the points defining the great circle
    pub fn end(&self) -> SphericalPoint {
        self.end
    }

    /// Returns the vector normal to the great circle
    pub fn normal(&self) -> Vector3<f32> {
        self.normal
    }

    /// Returns the intersections of this great circle and another one
    ///
    /// # Errors
    /// If the great circles are (essentially) parallel (equal to each other), returns [SphericalError::IdenticalGreatCircles] as then there is an infinite amount of intersections. You can handle this error as an equivalent of "all points on the circle are intersections".
    pub fn intersect_great_circle(&self, other: &GreatCircle) -> Result<[SphericalPoint; 2], SphericalError> {
        let normal1 = self.normal();
        let normal2 = other.normal();

        let res = normal1.cross(&normal2);
        if res.magnitude_squared() < VEC_LEN_IS_ZERO.powi(2) {
            return Err(SphericalError::IdenticalGreatCircles);
        }
        let res_norm = res.normalize();
        Ok([SphericalPoint::from_cartesian_vector3(res_norm), SphericalPoint::from_cartesian_vector3(-res_norm)])
    }

    /// Checks if the great circle contains the provided point
    pub fn contains_point(&self, point: &SphericalPoint) -> bool {
        let tolerance: f32 = 10e-6;
        let normal = self.start.cartesian().cross(&self.end.cartesian());
        normal.dot(&point.cartesian()) < tolerance
    }
}

#[cfg(feature = "serde")]
impl Serialize for GreatCircle {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let data = (self.start, self.end);
        data.serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for GreatCircle {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let (start, end) = <(SphericalPoint, SphericalPoint)>::deserialize(deserializer)?;
        GreatCircle::new(start, end).map_err(|e| serde::de::Error::custom(format!("invalid great circle: {:?}", e)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn gecaa_2020_theory_7() {
        let delta = 10e-2;

        let start_1 = SphericalPoint::new(-PI / 2.0, 0.0);
        let end_1 = SphericalPoint::new(0.0, 15.0 * PI / 180.0);
        let circle_1 = GreatCircle::new(start_1, end_1).expect("The points are fairly far away");

        let start_2 = SphericalPoint::new(-210.0 * PI / 180.0, 23.5 * PI / 180.0); // Switch RA direction as the question measures azimuth from north to east
        let end_2 = SphericalPoint::new(-255.0 * PI / 180.0, 75.0 * PI / 180.0); // Switch RA direction as the question measures azimuth from north to east
        let circle_2 = GreatCircle::new(start_2, end_2).expect("The points are fairly far away");

        let intersections = circle_1.intersect_great_circle(&circle_2).expect("The paths are not parallel");

        #[cfg(test)]
        dbg!(intersections);

        let [(ra_1, dec_1), (ra_2, dec_2)] = if intersections[0].ra() < intersections[1].ra() {
            [(intersections[1].ra(), intersections[1].dec()), (intersections[0].ra(), intersections[0].dec())]
        } else {
            [(intersections[0].ra(), intersections[0].dec()), (intersections[1].ra(), intersections[1].dec())]
        };

        let (ra_1_corr, dec_1_corr) = ((360.0 - 21.94) * PI / 180.0, 13.96 * PI / 180.0); // Once again switch RA direction as the question measures azimuth from north to east
        let (ra_2_corr, dec_2_corr) = ((360.0 - 201.94) * PI / 180.0, -13.96 * PI / 180.0); // Once again switch RA direction as the question measures azimuth from north to east

        assert!((ra_1_corr - ra_1).abs() < delta && (dec_1_corr - dec_1).abs() < delta);
        assert!((ra_2_corr - ra_2).abs() < delta && (dec_2_corr - dec_2).abs() < delta);
    }

    #[test]
    fn contains_point() {
        let equator = GreatCircle::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, 0.0)).expect("The points are far enough");
        let point_on_equator = SphericalPoint::new(PI / 3.0, 0.0);
        let point_outside_equator = SphericalPoint::new(0.0, 10e-5);

        assert!(equator.contains_point(&point_on_equator));
        assert!(!equator.contains_point(&point_outside_equator));

        let circle_2 = GreatCircle::new(SphericalPoint::new(PI / 5.0, PI / 4.0), SphericalPoint::new(PI / 5.0, -PI / 4.0)).expect("The points are far enough");
        assert!(circle_2.contains_point(&SphericalPoint::new(PI / 5.0, -PI / 7.0)));
    }

    #[test]
    fn perpendicular_through_point() {
        let tolerance = 10e-6;

        let equator = GreatCircle::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(-PI / 2.0, 0.0)).expect("The points are fairly far away");
        let greenwich_meridian = equator.perpendicular_through_point(&SphericalPoint::new(0.0, 0.0)).expect("The points are fairly far away");
        let new_pole = SphericalPoint::from_cartesian_vector3(greenwich_meridian.normal());
        let new_pole_corr_1 = SphericalPoint::new(PI / 2.0, 0.0);
        let new_pole_corr_2 = SphericalPoint::new(-PI / 2.0, 0.0);
        assert!(new_pole.approximately_equals(&new_pole_corr_1, tolerance) || new_pole.approximately_equals(&new_pole_corr_2, tolerance));

        let equator = GreatCircle::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(-PI / 2.0, 0.0)).expect("The points are fairly far away");
        let greenwich_meridian = equator.perpendicular_through_point(&SphericalPoint::new(0.0, PI / 2.0)); // Telling it to pass through the North Pole, which it should always do
        assert!(greenwich_meridian.is_err());

        let circle_1 = GreatCircle::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, PI / 6.0)).expect("The points are fairly far away");
        let circle_2 = circle_1.perpendicular_through_point(&SphericalPoint::new(0.0, PI / 2.0)).expect("The points are fairly far away");
        let new_pole = SphericalPoint::from_cartesian_vector3(circle_2.normal());
        let new_pole_corr_1 = SphericalPoint::new(0.0, 0.0);
        let new_pole_corr_2 = SphericalPoint::new(PI, 0.0);
        assert!(new_pole.approximately_equals(&new_pole_corr_1, tolerance) || new_pole.approximately_equals(&new_pole_corr_2, tolerance));

        let circle_1 = GreatCircle::new(SphericalPoint::new(0.0, 0.0), SphericalPoint::new(PI / 2.0, PI / 6.0)).expect("The points are fairly far away");
        let circle_2 = circle_1.perpendicular_through_point(&SphericalPoint::new(0.0, 0.0)).expect("The points are fairly far away");
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
            let orig = GreatCircle::new(point1, point2).expect("Valid points should create a GreatCircle");
            let ser = serde_json::to_string(&orig).expect("Serialization failed");
            let deser: GreatCircle = serde_json::from_str(&ser).expect("Deserialization failed");

            assert!((orig.start().ra() - deser.start().ra()).abs() < f32::EPSILON, "Start RA values do not match");
            assert!((orig.start().dec() - deser.start().dec()).abs() < f32::EPSILON, "Start Dec values do not match");
            assert!((orig.end().ra() - deser.end().ra()).abs() < f32::EPSILON, "End RA values do not match");
            assert!((orig.end().dec() - deser.end().dec()).abs() < f32::EPSILON, "End Dec values do not match");
            assert!((orig.normal() - deser.normal()).magnitude() < f32::EPSILON, "Normal vector does not match after deserialization");
        }
    }
}
