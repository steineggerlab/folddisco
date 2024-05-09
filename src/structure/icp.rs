// AI generated. IMPORTANT: TODO: NEED TO BE CHECKED.
use nalgebra::{DMatrix, DVector, U3, SVD};
use crate::structure::coordinate::Coordinate;

pub fn icp(source: &Vec<Coordinate>, target: &Vec<Coordinate>, max_iterations: usize) -> (DMatrix<f64>, DVector<f64>) {
    let mut T = DMatrix::<f64>::identity(3, 3);
    let mut t = DVector::<f64>::zeros(3);

    for _ in 0..max_iterations {
        let mut closest_points = Vec::new();

        for point in source {
            let closest_point = target.iter().min_by_key(|&r| r.distance_squared(point)).unwrap();
            closest_points.push(*closest_point);
        }

        let source_mean = mean(source);
        let target_mean = mean(&closest_points);

        let source_demean: Vec<_> = source.iter().map(|&p| p - source_mean).collect();
        let target_demean: Vec<_> = closest_points.iter().map(|&p| p - target_mean).collect();

        let H = source_demean.iter().zip(target_demean.iter()).map(|(&a, &b)| a.outer(b)).sum();

        let svd = SVD::new(H, true, true);
        let R = svd.u.unwrap() * svd.v.unwrap().transpose();
        let t = target_mean - R * source_mean;

        T = R * T;
        t = R * t + t;

        // Transform source points
        for point in source {
            *point = T * *point + t;
        }
    }

    (T, t)
}

fn mean(points: &Vec<Coordinate>) -> Coordinate {
    points.iter().sum::<Coordinate>() / (points.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn test_icp() {
        let source = vec![
            Coordinate::new(1.0, 2.0, 3.0),
            Coordinate::new(4.0, 5.0, 6.0),
            Coordinate::new(7.0, 8.0, 9.0),
        ];

        let target = vec![
            Coordinate::new(2.0, 4.0, 6.0),
            Coordinate::new(5.0, 7.0, 9.0),
            Coordinate::new(8.0, 10.0, 12.0),
        ];

        let (T, t) = icp(&source, &target, 10);

        // Check if the transformation matrix and translation vector are as expected.
        // These expected values are just for illustration and may not be correct.
        // You should replace them with the actual expected values.
        let expected_T = DMatrix::<f64>::identity(3, 3);
        let expected_t = DVector::<f64>::new(1.0, 2.0, 3.0);

        assert_eq!(T, expected_T);
        assert_eq!(t, expected_t);
    }
}