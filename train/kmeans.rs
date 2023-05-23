// Generated
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

fn main() {
    // Read data
    let file = File::open("yeast_raw_feature_100000.tsv").unwrap();
    let reader = BufReader::new(file);
    let mut data = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let mut row = Vec::new();
        for val in line.split("\t") {
            row.push(val.parse::<f32>().unwrap());
        }
        data.push(row);
    }

    // K means clustering
    let k = 1024;
    let mut centroids = Vec::new();
    for i in 0..k {
        centroids.push(data[i].clone());
    }
    loop {
        let mut clusters = vec![Vec::new(); k];
        for row in &data {
            let mut min_distance = std::f32::INFINITY;
            let mut closest_centroid = 0;
            for (i, centroid) in centroids.iter().enumerate() {
                let distance = euclidean_distance(row, centroid);
                if distance < min_distance {
                    min_distance = distance;
                    closest_centroid = i;
                }
            }
            clusters[closest_centroid].push(row.clone());
        }
        let mut new_centroids = Vec::new();
        for cluster in &clusters {
            let mut centroid = vec![0.0; 6];
            for row in cluster {
                for (i, val) in row.iter().enumerate() {
                    centroid[i] += val;
                }
            }
            for val in &mut centroid {
                *val /= cluster.len() as f32;
            }
            new_centroids.push(centroid);
        }
        if centroids == new_centroids {
            // Save clusters into file
            let mut file = File::create("yeast_clusters.tsv").unwrap();
            for (i, cluster) in clusters.iter().enumerate() {
                file.write(format!("{}\t", i).as_bytes()).unwrap();
                for row in cluster {
                    for val in row {
                        file.write(format!("{}\t", val).as_bytes()).unwrap();
                    }
                    file.write("\n".as_bytes()).unwrap();
                }
                file.write("\n".as_bytes()).unwrap();
            }
            break;

        }
        centroids = new_centroids;
    }

    // Save centroids into file
    let mut file = File::create("yeast_centroids.tsv").unwrap();
    for (i, centroid) in centroids.iter().enumerate() {
        file.write(format!("{}\t", i).as_bytes()).unwrap();
        for val in centroid {
            file.write(format!("{}\t", val).as_bytes()).unwrap();
        }
        file.write("\n".as_bytes()).unwrap();
    }

}

fn euclidean_distance(a: &[f32], b: &[f32]) -> f32 {
    let mut sum = 0.0;
    for i in 0..a.len() {
        sum += (a[i] - b[i]).powi(2);
    }
    sum.sqrt()
}

fn fast_assign_raw_val_to_cluster(raw_val: f32, centroids: &[Vec<f32>]) -> usize {
    let mut min_distance = std::f32::INFINITY;
    let mut closest_centroid = 0;
    for (i, centroid) in centroids.iter().enumerate() {
        let distance = (raw_val - centroid[0]).powi(2);
        if distance < min_distance {
            min_distance = distance;
            closest_centroid = i;
        }
    }
    closest_centroid
}