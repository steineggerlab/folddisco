// AI generated. IMPORTANT: TODO: NEED TO BE CHECKED.
use nalgebra::{DMatrix, DVector, U3, SVD};
use crate::structure::coordinate::Coordinate;

const DEFAULT_MAX_ITERATIONS: usize = 50;

pub fn icp(source: &Vec<Coordinate>, target: &Vec<Coordinate>, max_iterations: usize) -> (DMatrix<f32>, DVector<f32>) {
    let mut rotation = DMatrix::<f32>::identity(3, 3);
    let mut translation = DVector::<f32>::zeros(3);
    
    let mut transformed_source = source.clone();
    let mut prev_rmsd: f32 = std::f32::MAX;
    let rmsd_change_threshold: f32 = 0.001;
    for i in 0..max_iterations {
        let mut closest_points = Vec::new();

        for point in &transformed_source {
            let closest_point = target.iter().min_by(|&a, &b| {
                let distance_a = point.distance(a);
                let distance_b = point.distance(b);
                distance_a.partial_cmp(&distance_b).unwrap()
            });
            closest_points.push(*closest_point.unwrap());
        }

        let source_mean = mean(&transformed_source);
        let target_mean = mean(&closest_points);

        let source_demean: Vec<Coordinate> = transformed_source.iter().map(|&p| p.sub(&source_mean)).collect();
        let target_demean: Vec<Coordinate> = closest_points.iter().map(|&p| p.sub(&target_mean)).collect();
        // Covariance matrix
        let H = source_demean.iter().zip(target_demean.iter())
            .fold(DMatrix::<f32>::zeros(3, 3), |acc, (source_point, target_point)| {
                let source_matrix = coordinate_to_matrix(source_point);
                let target_matrix = coordinate_to_matrix(target_point);
                acc + source_matrix * target_matrix.transpose()
            });

        
        let svd = SVD::new(H, true, true);
        let mut u = svd.u.unwrap();
        let v_t = svd.v_t.unwrap();
        let mut R = &u * &v_t;
        if R.determinant() < 0.0 {
            u.column_mut(2).scale_mut(-1.0);
            R = &u * &v_t;
        }
        let target_mean_vec = coordinate_to_matrix(&target_mean);
        let source_mean_vec = coordinate_to_matrix(&source_mean);
        let t = target_mean_vec - &R * source_mean_vec;

        // Update transformation matrix and transform source points
        rotation = &R * rotation;
        translation = &R * &translation + &t;
        transformed_source = transform_coordinate_vector(&rotation, &translation, &source);

        // Calculate RMSD
        let mut rmsd = transformed_source.iter().zip(target.iter())
            .map(|(&transformed_point, &target_point)| (transformed_point.sub(&target_point)).norm().powi(2))
            .sum::<f32>() / (transformed_source.len() as f32);
        rmsd = rmsd.sqrt();
        println!("RMSD: {}", rmsd);
        if (prev_rmsd - rmsd).abs() < rmsd_change_threshold {
            break;
        }
        prev_rmsd = rmsd;

        println!("Transformed source: {:?} at iteration {}; cost={}", transformed_source, i, rmsd);
    }

    (rotation, translation)
}
fn mean(points: &Vec<Coordinate>) -> Coordinate {
    points.iter().fold(Coordinate::new(0.0, 0.0, 0.0), |acc, &p| acc.add(&p)).scale(1.0 / points.len() as f32)
}

fn coordinate_to_matrix(point: &Coordinate) -> DMatrix<f32> {
    DMatrix::<f32>::from_vec(3, 1, vec![point.x, point.y, point.z])
}
fn coordinate_vector_to_matrix(points: &Vec<Coordinate>) -> DMatrix<f32> {
    let mut matrix = DMatrix::<f32>::zeros(3, points.len());
    for (i, point) in points.iter().enumerate() {
        matrix[(0, i)] = point.x;
        matrix[(1, i)] = point.y;
        matrix[(2, i)] = point.z;
    }
    matrix
}

fn matrix_to_coordinate(matrix: &DMatrix<f32>) -> Coordinate {
    Coordinate::new(matrix[(0, 0)], matrix[(1, 0)], matrix[(2, 0)])
}

fn transform_coordinate_vector(rotation: &DMatrix<f32>, translation: &DVector<f32>, points: &Vec<Coordinate>) -> Vec<Coordinate> {
    let points_matrix = coordinate_vector_to_matrix(points);
    let transformed_points_matrix = rotation * points_matrix;
    let mut transformed_points = Vec::new();
    for i in 0..points.len() {
        let mut transformed_point = Coordinate::new(
            transformed_points_matrix[(0, i)],
            transformed_points_matrix[(1, i)],
            transformed_points_matrix[(2, i)],
        );
        transformed_point.x += translation[0];
        transformed_point.y += translation[1];
        transformed_point.z += translation[2];
        transformed_points.push(transformed_point);
    }
    transformed_points
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

        
        for i in 1..=10 {
            let (T, t) = icp(&source, &target, i);
            println!("T: {:?}", T);
            println!("t: {:?}", t);
            let transformed_points = transform_coordinate_vector(&T, &t, &source);
            for point in &transformed_points {
                println!("{:?}", point);
            }
            let rmsd = transformed_points.iter().zip(target.iter()).fold(
                0.0, 
                |acc, (p1, p2)| {
                    let squared_distance = p1.distance(p2).powi(2);
                    acc + squared_distance / transformed_points.len() as f32
                }
            ).sqrt();
            println!("RMSD: {}", rmsd);
        }
    }
    
    #[test]
    fn test_actual_coordinate() {
        let source = vec![
            Coordinate::new(6.994, 8.354, 42.405),
            Coordinate::new(9.429, 7.479, 48.266),
            Coordinate::new(5.547, 0.158, 42.050),
        ];
        let target1 = vec![
            Coordinate::new(-12.833, 3.134, -7.780),
            Coordinate::new(-5.720, -2.218, -3.368),
            Coordinate::new(-13.958, -1.741, -4.223),
        ];
        let target2 = vec![
            Coordinate::new(-13.958, -1.741, -4.223),
            Coordinate::new(-12.833, 3.134, -7.780),
            Coordinate::new(-5.720, -2.218, -3.368),
        ];
        let iter_vec: Vec<usize> = vec![50];
        for i in iter_vec {
            let (T1, t1) = icp(&source, &target1, i);
            let transformed_points1 = transform_coordinate_vector(&T1, &t1, &source);
            let rmsd1 = transformed_points1.iter().zip(target1.iter()).fold(
                0.0, 
                |acc, (p1, p2)| {
                    let squared_distance = p1.distance(p2).powi(2);
                    acc + squared_distance / transformed_points1.len() as f32
                }
            ).sqrt();
            println!("Iteration: {}, Transformed Points: {:?}, RMSD1: {}", i, transformed_points1, rmsd1);

            let (T2, t2) = icp(&source, &target2, i);
            let transformed_points2 = transform_coordinate_vector(&T2, &t2, &source);
            let rmsd2 = transformed_points2.iter().zip(target2.iter()).fold(
                0.0, 
                |acc, (p1, p2)| {
                    let squared_distance = p1.distance(p2).powi(2);
                    acc + squared_distance / transformed_points2.len() as f32
                }
            ).sqrt();
            println!("Iteration: {}, Transformed Points: {:?}, RMSD2: {}", i, transformed_points2, rmsd2);
        }
    }
    
}