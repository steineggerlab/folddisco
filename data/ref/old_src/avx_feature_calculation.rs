
pub fn compute_pairwise_features_old(compact: &CompactStructure) -> Vec<[f32; 6]> {
    let n = compact.num_residues;
    let mut features = Vec::with_capacity(n * n); // Reserve space for N * N features

    unsafe {
        for i in 0..n {
            let ca1 = &compact.ca_vector2;
            let cb1 = &compact.cb_vector2;
            let n1 = &compact.n_vector2;

            // Load the first vectors for ca1, cb1, and n1
            let ca1x = _mm256_set1_ps(ca1.0[i]);
            let ca1y = _mm256_set1_ps(ca1.1[i]);
            let ca1z = _mm256_set1_ps(ca1.2[i]);
            let cb1x = _mm256_set1_ps(cb1.0[i]);
            let cb1y = _mm256_set1_ps(cb1.1[i]);
            let cb1z = _mm256_set1_ps(cb1.2[i]);
            let n1x = _mm256_set1_ps(n1.0[i]);
            let n1y = _mm256_set1_ps(n1.1[i]);
            let n1z = _mm256_set1_ps(n1.2[i]);

            for j in (0..n).step_by(8) {
                // Load coordinates for 8 pairs of residues
                let ca2x = _mm256_loadu_ps(&ca1.0[j]);
                let ca2y = _mm256_loadu_ps(&ca1.1[j]);
                let ca2z = _mm256_loadu_ps(&ca1.2[j]);
                // let cb2x = _mm256_loadu_ps(&cb1.0[j]);
                // let cb2y = _mm256_loadu_ps(&cb1.1[j]);
                // let cb2z = _mm256_loadu_ps(&cb1.2[j]);

                // Compute squared differences and distances for CA
                let cadx = _mm256_sub_ps(ca2x, ca1x);
                let cady = _mm256_sub_ps(ca2y, ca1y);
                let cadz = _mm256_sub_ps(ca2z, ca1z);

                let ca2dx = _mm256_mul_ps(cadx, cadx);
                let ca2dy = _mm256_mul_ps(cady, cady);
                let ca2dz = _mm256_mul_ps(cadz, cadz);

                let cadist2 = _mm256_add_ps(ca2dx, _mm256_add_ps(ca2dy, ca2dz));
                // let ca_dist = _mm256_sqrt_ps(cadist2);

                // Skip elements where ca_dist > 20.0
                let mask = _mm256_cmp_ps(cadist2, _mm256_set1_ps(400.0), _CMP_LE_OS); // Compare ca_dist <= 20.0
                if _mm256_testz_ps(mask, mask) != 0 {
                    continue; // Skip iteration if all distances are greater than 20.0
                }
                let ca_dist = _mm256_sqrt_ps(cadist2);

                // Load CB and N2 coordinates if ca_dist <= 20.0
                let cb2x = _mm256_loadu_ps(&cb1.0[j]);
                let cb2y = _mm256_loadu_ps(&cb1.1[j]);
                let cb2z = _mm256_loadu_ps(&cb1.2[j]);
                let n2x = _mm256_loadu_ps(&n1.0[j]);
                let n2y = _mm256_loadu_ps(&n1.1[j]);
                let n2z = _mm256_loadu_ps(&n1.2[j]);

                // Compute squared differences and distances for CB
                let cbdx = _mm256_sub_ps(cb2x, cb1x);
                let cbdy = _mm256_sub_ps(cb2y, cb1y);
                let cbdz = _mm256_sub_ps(cb2z, cb1z);

                let cb2dx = _mm256_mul_ps(cbdx, cbdx);
                let cb2dy = _mm256_mul_ps(cbdy, cbdy);
                let cb2dz = _mm256_mul_ps(cbdz, cbdz);

                let cbdist2 = _mm256_add_ps(cb2dx, _mm256_add_ps(cb2dy, cb2dz));
                let cb_dist = _mm256_sqrt_ps(cbdist2);

                // Continue with angle and torsion computations
                // Compute angle between ca1-cb1 vector and ca2-cb2 vector
                let ca_cb1x = _mm256_sub_ps(cb1x, ca1x);
                let ca_cb1y = _mm256_sub_ps(cb1y, ca1y);
                let ca_cb1z = _mm256_sub_ps(cb1z, ca1z);
                let ca_cb2x = _mm256_sub_ps(cb2x, ca2x);
                let ca_cb2y = _mm256_sub_ps(cb2y, ca2y);
                let ca_cb2z = _mm256_sub_ps(cb2z, ca2z);

                let angle_cos = _mm256_div_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(ca_cb1x, ca_cb2x),
                        _mm256_add_ps(
                            _mm256_mul_ps(ca_cb1y, ca_cb2y),
                            _mm256_mul_ps(ca_cb1z, ca_cb2z)
                        )
                    ),
                    _mm256_mul_ps(
                        _mm256_sqrt_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(ca_cb1x, ca_cb1x),
                                    _mm256_mul_ps(ca_cb1y, ca_cb1y)
                                ),
                                _mm256_mul_ps(ca_cb1z, ca_cb1z)
                            )
                        ),
                        _mm256_sqrt_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(ca_cb2x, ca_cb2x),
                                    _mm256_mul_ps(ca_cb2y, ca_cb2y)
                                ),
                                _mm256_mul_ps(ca_cb2z, ca_cb2z)
                            )
                        )
                    )
                );

                let angle_sin = _mm256_sqrt_ps(_mm256_sub_ps(_mm256_set1_ps(1.0), _mm256_mul_ps(angle_cos, angle_cos)));

                // Continue with torsion angle computations (similar to existing code)
                let torsion_angle1_tan = compute_torsion_angle(cb1x, cb1y, cb1z, cb2x, cb2y, cb2z, ca1x, ca1y, ca1z, n1x, n1y, n1z);
                let torsion_angle2_tan = compute_torsion_angle(cb1x, cb1y, cb1z, cb2x, cb2y, cb2z, ca2x, ca2y, ca2z, n2x, n2y, n2z);

                // Pack results and store
                let mut temp_ca_dist = [0.0; 8];
                let mut temp_cb_dist = [0.0; 8];
                let mut temp_angle_cos = [0.0; 8];
                let mut temp_angle_sin = [0.0; 8];
                let mut temp_torsion1 = [0.0; 8];
                let mut temp_torsion2 = [0.0; 8];

                _mm256_storeu_ps(temp_ca_dist.as_mut_ptr(), ca_dist);
                _mm256_storeu_ps(temp_cb_dist.as_mut_ptr(), cb_dist);
                _mm256_storeu_ps(temp_angle_cos.as_mut_ptr(), angle_cos);
                _mm256_storeu_ps(temp_angle_sin.as_mut_ptr(), angle_sin);
                _mm256_storeu_ps(temp_torsion1.as_mut_ptr(), torsion_angle1_tan);
                _mm256_storeu_ps(temp_torsion2.as_mut_ptr(), torsion_angle2_tan);

                for k in 0..8 {
                    if temp_ca_dist[k] <= 20.0 {
                        features.push([temp_ca_dist[k], temp_cb_dist[k], temp_angle_cos[k], temp_angle_sin[k], temp_torsion1[k].atan(), temp_torsion2[k].atan()]);
                    }
                }
            }
        }
    }
    features
}

    
// Function to compute pairwise features with SIMD
// TEMP
fn compute_pairwise_features_temp(compact: &CompactStructure) -> Vec<[f32; 6]> {
    // SIMD acceleration for pdbtr pairwise feature computation. 8 x f32 with avx1.0
    // Return vector of ca_dist, cb_dist, angle_cos, and angle_sin for now
    let n = compact.num_residues;
    let mut features = Vec::with_capacity(n * n); // Reserve space for N * N features

    unsafe {
        for i in 0..n {
            let ca1 = &compact.ca_vector2;
            let cb1 = &compact.cb_vector2;
            let n1 = &compact.n_vector2;
            let ca1x = _mm256_set1_ps(ca1.0[i]);
            let ca1y = _mm256_set1_ps(ca1.1[i]);
            let ca1z = _mm256_set1_ps(ca1.2[i]);
            let cb1x = _mm256_set1_ps(cb1.0[i]);
            let cb1y = _mm256_set1_ps(cb1.1[i]);
            let cb1z = _mm256_set1_ps(cb1.2[i]);
            let n1x = _mm256_set1_ps(n1.0[i]);
            let n1y = _mm256_set1_ps(n1.1[i]);
            let n1z = _mm256_set1_ps(n1.2[i]);
            
            for j in (0..n).step_by(8) {

                // Load x2, y2, z2 for 8 pairs of coordinates
                let ca2x = _mm256_loadu_ps(&ca1.0[j]);
                let ca2y = _mm256_loadu_ps(&ca1.1[j]);
                let ca2z = _mm256_loadu_ps(&ca1.2[j]);
                let cb2x = _mm256_loadu_ps(&cb1.0[j]);
                let cb2y = _mm256_loadu_ps(&cb1.1[j]);
                let cb2z = _mm256_loadu_ps(&cb1.2[j]);
                let n2x = _mm256_loadu_ps(&n1.0[j]);
                let n2y = _mm256_loadu_ps(&n1.1[j]);
                let n2z = _mm256_loadu_ps(&n1.2[j]);
                
                // Compute the squared differences for distance calculation
                let cadx = _mm256_sub_ps(ca2x, ca1x);
                let cady = _mm256_sub_ps(ca2y, ca1y);
                let cadz = _mm256_sub_ps(ca2z, ca1z);
                
                let ca2dx = _mm256_mul_ps(cadx, cadx);
                let ca2dy = _mm256_mul_ps(cady, cady);
                let ca2dz = _mm256_mul_ps(cadz, cadz);

                let cbdx = _mm256_sub_ps(cb2x, cb1x);
                let cbdy = _mm256_sub_ps(cb2y, cb1y);
                let cbdz = _mm256_sub_ps(cb2z, cb1z);
                
                let cb2dx = _mm256_mul_ps(cbdx, cbdx);
                let cb2dy = _mm256_mul_ps(cbdy, cbdy);
                let cb2dz = _mm256_mul_ps(cbdz, cbdz);
                
                let cadist2 = _mm256_add_ps(ca2dx, _mm256_add_ps(ca2dy, ca2dz));
                let cbdist2 = _mm256_add_ps(cb2dx, _mm256_add_ps(cb2dy, cb2dz));
                let ca_dist = _mm256_sqrt_ps(cadist2);
                let cb_dist = _mm256_sqrt_ps(cbdist2);
                
                // Compute angle between ca1-cb1 vector and ca2-cb2 vector
                let ca_cb1x = _mm256_sub_ps(cb1x, ca1x);
                let ca_cb1y = _mm256_sub_ps(cb1y, ca1y);
                let ca_cb1z = _mm256_sub_ps(cb1z, ca1z);
                let ca_cb2x = _mm256_sub_ps(cb2x, ca2x);
                let ca_cb2y = _mm256_sub_ps(cb2y, ca2y);
                let ca_cb2z = _mm256_sub_ps(cb2z, ca2z);
                
                let angle_cos = _mm256_div_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(ca_cb1x, ca_cb2x),
                        _mm256_add_ps(
                            _mm256_mul_ps(ca_cb1y, ca_cb2y),
                            _mm256_mul_ps(ca_cb1z, ca_cb2z)
                        )
                    ),
                    _mm256_mul_ps(
                        _mm256_sqrt_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(ca_cb1x, ca_cb1x),
                                    _mm256_mul_ps(ca_cb1y, ca_cb1y)
                                ),
                                _mm256_mul_ps(ca_cb1z, ca_cb1z)
                            )
                        ),
                        _mm256_sqrt_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(ca_cb2x, ca_cb2x),
                                    _mm256_mul_ps(ca_cb2y, ca_cb2y)
                                ),
                                _mm256_mul_ps(ca_cb2z, ca_cb2z)
                            )
                        )
                    )
                );
                let angle_sin = _mm256_sqrt_ps(_mm256_sub_ps(_mm256_set1_ps(1.0), _mm256_mul_ps(angle_cos, angle_cos)));
                
                // // TODO: Calculate the torsion angles
                // // First torsion angle: n1-ca1-cb1-cb2
                // // Second torsion angle: cb1-cb2-ca2-n2
                let n1ca1x = _mm256_sub_ps(ca1x, n1x);
                let n1ca1y = _mm256_sub_ps(ca1y, n1y);
                let n1ca1z = _mm256_sub_ps(ca1z, n1z);
                let n2ca2x = _mm256_sub_ps(ca2x, n2x);
                let n2ca2y = _mm256_sub_ps(ca2y, n2y);
                let n2ca2z = _mm256_sub_ps(ca2z, n2z);
                
                // FILL IN THIS PART
                                // First torsion angle: n1-ca1-cb1-cb2
                let v1x = _mm256_sub_ps(cb1x, ca1x);
                let v1y = _mm256_sub_ps(cb1y, ca1y);
                let v1z = _mm256_sub_ps(cb1z, ca1z);

                let v2x = _mm256_sub_ps(ca1x, n1x);
                let v2y = _mm256_sub_ps(ca1y, n1y);
                let v2z = _mm256_sub_ps(ca1z, n1z);

                let v3x = _mm256_sub_ps(cb2x, cb1x);
                let v3y = _mm256_sub_ps(cb2y, cb1y);
                let v3z = _mm256_sub_ps(cb2z, cb1z);

                // Cross product v1 x v2
                let cross1x = _mm256_sub_ps(_mm256_mul_ps(v1y, v2z), _mm256_mul_ps(v1z, v2y));
                let cross1y = _mm256_sub_ps(_mm256_mul_ps(v1z, v2x), _mm256_mul_ps(v1x, v2z));
                let cross1z = _mm256_sub_ps(_mm256_mul_ps(v1x, v2y), _mm256_mul_ps(v1y, v2x));

                // Cross product v2 x v3
                let cross2x = _mm256_sub_ps(_mm256_mul_ps(v2y, cbdz), _mm256_mul_ps(v2z, cbdy));
                let cross2y = _mm256_sub_ps(_mm256_mul_ps(v2z, cbdx), _mm256_mul_ps(v2x, cbdz));
                let cross2z = _mm256_sub_ps(_mm256_mul_ps(v2x, cbdy), _mm256_mul_ps(v2y, cbdx));

                // Dot product of cross products
                let dot_cross = _mm256_add_ps(
                    _mm256_mul_ps(cross1x, cross2x),
                    _mm256_add_ps(
                        _mm256_mul_ps(cross1y, cross2y),
                        _mm256_mul_ps(cross1z, cross2z)
                    )
                );

                // Compute magnitudes
                let cross1_mag = _mm256_sqrt_ps(_mm256_add_ps(
                    _mm256_mul_ps(cross1x, cross1x),
                    _mm256_add_ps(_mm256_mul_ps(cross1y, cross1y), _mm256_mul_ps(cross1z, cross1z))
                ));
                let cross2_mag = _mm256_sqrt_ps(_mm256_add_ps(
                    _mm256_mul_ps(cross2x, cross2x),
                    _mm256_add_ps(_mm256_mul_ps(cross2y, cross2y), _mm256_mul_ps(cross2z, cross2z))
                ));

                // Calculate torsion angle. cannot use atan2. just return as tangent
                let torsion_angle1_tan = _mm256_div_ps(cross1_mag, dot_cross);
                // Second torsion angle: cb1-cb2-ca2-n2
                let v4x = _mm256_sub_ps(ca2x, cb2x);
                let v4y = _mm256_sub_ps(ca2y, cb2y);
                let v4z = _mm256_sub_ps(ca2z, cb2z);

                let v5x = _mm256_sub_ps(n2x, ca2x);
                let v5y = _mm256_sub_ps(n2y, ca2y);
                let v5z = _mm256_sub_ps(n2z, ca2z);

                // Cross product cbd x v4
                let cross3x = _mm256_sub_ps(_mm256_mul_ps(cbdy, v4z), _mm256_mul_ps(cbdz, v4y));
                let cross3y = _mm256_sub_ps(_mm256_mul_ps(cbdz, v4x), _mm256_mul_ps(cbdx, v4z));
                let cross3z = _mm256_sub_ps(_mm256_mul_ps(cbdx, v4y), _mm256_mul_ps(cbdy, v4x));
                
                // Cross product v4 x v5
                let cross4x = _mm256_sub_ps(_mm256_mul_ps(v4y, v5z), _mm256_mul_ps(v4z, v5y));
                let cross4y = _mm256_sub_ps(_mm256_mul_ps(v4z, v5x), _mm256_mul_ps(v4x, v5z));
                let cross4z = _mm256_sub_ps(_mm256_mul_ps(v4x, v5y), _mm256_mul_ps(v4y, v5x));
                
                // Dot product of cross products
                let dot_cross2 = _mm256_add_ps(
                    _mm256_mul_ps(cross3x, cross4x),
                    _mm256_add_ps(
                        _mm256_mul_ps(cross3y, cross4y),
                        _mm256_mul_ps(cross3z, cross4z)
                    )
                );
                // Compute magnitudes
                let cross3_mag = _mm256_sqrt_ps(_mm256_add_ps(
                    _mm256_mul_ps(cross3x, cross3x),
                    _mm256_add_ps(_mm256_mul_ps(cross3y, cross3y), _mm256_mul_ps(cross3z, cross3z))
                ));
                let cross4_mag = _mm256_sqrt_ps(_mm256_add_ps(
                    _mm256_mul_ps(cross4x, cross4x),
                    _mm256_add_ps(_mm256_mul_ps(cross4y, cross4y), _mm256_mul_ps(cross4z, cross4z))
                ));
                // Calculate torsion angle. cannot use atan2. just return as tangent
                let torsion_angle2_tan = _mm256_div_ps(cross3_mag, dot_cross2);
                // // Pack the integers into u32 and apply bit shifting if needed
                // let shifted_distances = _mm_slli_epi32(distance_int, 8); // Example bit shift

                // // Store the results
                // let mut temp_results = [0u32; 4];
                // _mm_storeu_si128(temp_results.as_mut_ptr() as *mut __m128i, shifted_distances);
                
                // Temporary lines for debugging: set values to 0
                // let ca_dist = _mm256_set1_ps(0.0);
                // let cb_dist = _mm256_set1_ps(0.0);
                // let angle_cos = _mm256_set1_ps(0.0);
                // let angle_sin = _mm256_set1_ps(0.0);
                // let torsion_angle1_tan = _mm256_set1_ps(0.0);
                // let torsion_angle2_tan = _mm256_set1_ps(0.0);
                
                
                
                
                let mut temp_ca_dist = [0.0; 8];
                let mut temp_cb_dist = [0.0; 8];
                let mut temp_angle_cos = [0.0; 8];
                let mut temp_angle_sin = [0.0; 8];
                let mut temp_torsion1 = [0.0; 8];
                let mut temp_torsion2 = [0.0; 8];
                _mm256_storeu_ps(temp_ca_dist.as_mut_ptr(), ca_dist);
                _mm256_storeu_ps(temp_cb_dist.as_mut_ptr(), cb_dist);
                _mm256_storeu_ps(temp_angle_cos.as_mut_ptr(), angle_cos);
                _mm256_storeu_ps(temp_angle_sin.as_mut_ptr(), angle_sin);
                _mm256_storeu_ps(temp_torsion1.as_mut_ptr(), torsion_angle1_tan);
                _mm256_storeu_ps(temp_torsion2.as_mut_ptr(), torsion_angle2_tan);
                for k in 0..8 {
                    features.push([temp_ca_dist[k], temp_cb_dist[k], temp_angle_cos[k], temp_angle_sin[k], temp_torsion1[k].atan(), temp_torsion2[k].atan()]);
                }
            }
        }
    }
    features
}

pub fn compute_pairwise_features(compact: &CompactStructure) -> Vec<u32> {
    let n = compact.num_residues;
    let mut features = Vec::with_capacity(n * n); // Reserve space for N * N features
    // Store results
    let mut temp_ca_dist = [0.0; 8];
    let mut temp_cb_dist = [0.0; 8];
    let mut temp_angle_cos = [0.0; 8];
    let mut temp_angle_sin = [0.0; 8];
    let mut temp_torsion1 = [0.0; 8];
    let mut temp_torsion2 = [0.0; 8];
    unsafe {
        for i in 0..n {
            let res1 = compact.get_res_name(i);
            let res1_u32 = map_aa_to_u8(res1) as u32;
            let ca1 = &compact.ca_vector2;
            let cb1 = &compact.cb_vector2;
            let n1 = &compact.n_vector2;

            // Load the first vectors for ca1, cb1, and n1
            let ca1x = _mm256_set1_ps(ca1.0[i]);
            let ca1y = _mm256_set1_ps(ca1.1[i]);
            let ca1z = _mm256_set1_ps(ca1.2[i]);
            let cb1x = _mm256_set1_ps(cb1.0[i]);
            let cb1y = _mm256_set1_ps(cb1.1[i]);
            let cb1z = _mm256_set1_ps(cb1.2[i]);
            let n1x = _mm256_set1_ps(n1.0[i]);
            let n1y = _mm256_set1_ps(n1.1[i]);
            let n1z = _mm256_set1_ps(n1.2[i]);
            
            for j in (0..n).step_by(8) {
            
                // Load coordinates for 8 pairs of residues
                let ca2x = _mm256_loadu_ps(&ca1.0[j]);
                let ca2y = _mm256_loadu_ps(&ca1.1[j]);
                let ca2z = _mm256_loadu_ps(&ca1.2[j]);

                // Compute squared differences and distances for CA and CB
                let ca_dist = compute_distance(ca1x, ca1y, ca1z, ca2x, ca2y, ca2z);
                // Skip elements where ca_dist > 20.0
                let mask = _mm256_cmp_ps(ca_dist, _mm256_set1_ps(20.0), _CMP_LE_OS); // Compare ca_dist <= 20.0
                if _mm256_testz_ps(mask, mask) != 0 {
                    continue; // Skip iteration if all distances are greater than 20.0
                }

                let cb2x = _mm256_loadu_ps(&cb1.0[j]);
                let cb2y = _mm256_loadu_ps(&cb1.1[j]);
                let cb2z = _mm256_loadu_ps(&cb1.2[j]);
                let n2x = _mm256_loadu_ps(&n1.0[j]);
                let n2y = _mm256_loadu_ps(&n1.1[j]);
                let n2z = _mm256_loadu_ps(&n1.2[j]);

                let res2_array = [
                    map_aa_to_u8(compact.get_res_name(j + 0)) as u32,
                    map_aa_to_u8(compact.get_res_name(j + 1)) as u32,
                    map_aa_to_u8(compact.get_res_name(j + 2)) as u32,
                    map_aa_to_u8(compact.get_res_name(j + 3)) as u32,
                    map_aa_to_u8(compact.get_res_name(j + 4)) as u32,
                    map_aa_to_u8(compact.get_res_name(j + 5)) as u32,
                    map_aa_to_u8(compact.get_res_name(j + 6)) as u32,
                    map_aa_to_u8(compact.get_res_name(j + 7)) as u32,
                ];
                
                let cb_dist = compute_distance(cb1x, cb1y, cb1z, cb2x, cb2y, cb2z);
                // Compute cos/sin angles for ca-cb vectors
                let (angle_cos, angle_sin) = compute_angle_cos_sin(ca1x, ca1y, ca1z, cb1x, cb1y, cb1z, ca2x, ca2y, ca2z, cb2x, cb2y, cb2z);
                
                // First torsion angle: n1-ca1-cb1-cb2
                let torsion_angle1_tan = compute_torsion_angle(ca1x, ca1y, ca1z, cb1x, cb1y, cb1z, cb2x, cb2y, cb2z, n1x, n1y, n1z);

                // Second torsion angle: cb1-cb2-ca2-n2
                let torsion_angle2_tan = compute_torsion_angle(cb1x, cb1y, cb1z, cb2x, cb2y, cb2z, ca2x, ca2y, ca2z, n2x, n2y, n2z);
                

                _mm256_storeu_ps(temp_ca_dist.as_mut_ptr(), ca_dist);
                _mm256_storeu_ps(temp_cb_dist.as_mut_ptr(), cb_dist);
                _mm256_storeu_ps(temp_angle_cos.as_mut_ptr(), angle_cos);
                _mm256_storeu_ps(temp_angle_sin.as_mut_ptr(), angle_sin);
                _mm256_storeu_ps(temp_torsion1.as_mut_ptr(), torsion_angle1_tan);
                _mm256_storeu_ps(temp_torsion2.as_mut_ptr(), torsion_angle2_tan);

                for k in 0..8 {
                    if temp_ca_dist[k] <= 20.0 {
                        let hash_value = pdb_tr::HashValue::perfect_hash_new2(res1_u32, res2_array[k], temp_ca_dist[k], temp_cb_dist[k], temp_angle_cos[k], temp_angle_sin[k], temp_torsion1[k], temp_torsion2[k]);
                        features.push(hash_value);
                    }
                }
            }
        }
    }
    features
}


// Helper to compute squared distances
#[inline(always)]
unsafe fn compute_distance(x1: __m256, y1: __m256, z1: __m256, x2: __m256, y2: __m256, z2: __m256) -> __m256 {
    let dx = _mm256_sub_ps(x2, x1);
    let dy = _mm256_sub_ps(y2, y1);
    let dz = _mm256_sub_ps(z2, z1);
    _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(dx, dx), _mm256_add_ps(_mm256_mul_ps(dy, dy), _mm256_mul_ps(dz, dz))))
}

// Helper to compute cos and sin of angles
#[inline(always)]
unsafe fn compute_angle_cos_sin(
    ca1x: __m256, ca1y: __m256, ca1z: __m256, cb1x: __m256, cb1y: __m256, cb1z: __m256,
    ca2x: __m256, ca2y: __m256, ca2z: __m256, cb2x: __m256, cb2y: __m256, cb2z: __m256
) -> (__m256, __m256) {
    let v1x = _mm256_sub_ps(cb1x, ca1x);
    let v1y = _mm256_sub_ps(cb1y, ca1y);
    let v1z = _mm256_sub_ps(cb1z, ca1z);
    let v2x = _mm256_sub_ps(cb2x, ca2x);
    let v2y = _mm256_sub_ps(cb2y, ca2y);
    let v2z = _mm256_sub_ps(cb2z, ca2z);

    let dot_product = _mm256_add_ps(_mm256_mul_ps(v1x, v2x), _mm256_add_ps(_mm256_mul_ps(v1y, v2y), _mm256_mul_ps(v1z, v2z)));

    let mag_v1 = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(v1x, v1x), _mm256_add_ps(_mm256_mul_ps(v1y, v1y), _mm256_mul_ps(v1z, v1z))));
    let mag_v2 = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(v2x, v2x), _mm256_add_ps(_mm256_mul_ps(v2y, v2y), _mm256_mul_ps(v2z, v2z))));

    let angle_cos = _mm256_div_ps(dot_product, _mm256_mul_ps(mag_v1, mag_v2));
    let angle_sin = _mm256_sqrt_ps(_mm256_sub_ps(_mm256_set1_ps(1.0), _mm256_mul_ps(angle_cos, angle_cos)));
    (angle_cos, angle_sin)
}


#[inline(always)]
unsafe fn compute_torsion_angle(
    x1: __m256, y1: __m256, z1: __m256, x2: __m256, y2: __m256, z2: __m256, 
    x3: __m256, y3: __m256, z3: __m256, x4: __m256, y4: __m256, z4: __m256
) -> __m256 {
    // Step 1: Compute v1, v2, and v3 (delta vectors)
    let v1x = _mm256_sub_ps(x2, x1);
    let v1y = _mm256_sub_ps(y2, y1);
    let v1z = _mm256_sub_ps(z2, z1);

    let v2x = _mm256_sub_ps(x3, x2);
    let v2y = _mm256_sub_ps(y3, y2);
    let v2z = _mm256_sub_ps(z3, z2);

    let v3x = _mm256_sub_ps(x4, x3);
    let v3y = _mm256_sub_ps(y4, y3);
    let v3z = _mm256_sub_ps(z4, z3);

    // Step 2: Compute cross products for v1 x v2 and v2 x v3
    let cross1x = _mm256_fmsub_ps(v1y, v2z, _mm256_mul_ps(v1z, v2y));
    let cross1y = _mm256_fmsub_ps(v1z, v2x, _mm256_mul_ps(v1x, v2z));
    let cross1z = _mm256_fmsub_ps(v1x, v2y, _mm256_mul_ps(v1y, v2x));

    let cross2x = _mm256_fmsub_ps(v2y, v3z, _mm256_mul_ps(v2z, v3y));
    let cross2y = _mm256_fmsub_ps(v2z, v3x, _mm256_mul_ps(v2x, v3z));
    let cross2z = _mm256_fmsub_ps(v2x, v3y, _mm256_mul_ps(v2y, v3x));

    // Step 3: Compute dot product of cross1 and cross2
    let dot_cross = _mm256_fmadd_ps(cross1x, cross2x, _mm256_fmadd_ps(cross1y, cross2y, _mm256_mul_ps(cross1z, cross2z)));

    // Step 4: Compute squared magnitudes of cross1 and cross2 (use the squared values to avoid recalculating later)
    let cross1_mag_sq = _mm256_fmadd_ps(cross1x, cross1x, _mm256_fmadd_ps(cross1y, cross1y, _mm256_mul_ps(cross1z, cross1z)));
    let cross2_mag_sq = _mm256_fmadd_ps(cross2x, cross2x, _mm256_fmadd_ps(cross2y, cross2y, _mm256_mul_ps(cross2z, cross2z)));

    // Step 5: Compute magnitudes (square roots)
    let mag_cross1 = _mm256_sqrt_ps(cross1_mag_sq);
    let mag_cross2 = _mm256_sqrt_ps(cross2_mag_sq);

    // Step 6: Multiply the magnitudes of cross1 and cross2
    let mag_product = _mm256_mul_ps(mag_cross1, mag_cross2);

    // Step 7: Compute the tangent of the torsion angle (magnitude of cross products divided by dot product)
    _mm256_div_ps(mag_product, dot_cross)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::core::CompactStructure;
    use crate::geometry::core::HashType;
    use crate::utils::combination::CombinationIterator;
    use crate::prelude::PDBReader;

    #[test]
    fn test_simd_feature() {
        let start = std::time::Instant::now();
        let pdb = "data/serine_peptidases_filtered/4cha.pdb";
        let pdb_reader = PDBReader::from_file(pdb).unwrap();
        let compact = pdb_reader.read_structure().unwrap().to_compact();
        let duration = start.elapsed();
        println!("Reading structure: {:?}", duration);
        println!("Number of residues: {}", compact.num_residues);
        println!("Square of number of residues: {}", compact.num_residues * compact.num_residues);
        let start = std::time::Instant::now();
        for i in 0..1 {
            let features = compute_pairwise_features(&compact);
        }
        let duration = start.elapsed();
        println!("SIMD feature computation: {:?}", duration);
        let features = compute_pairwise_features(&compact);
        // Print first 10 features
        println!("Number of features: {}", features.len());
        for i in 0..10 {
            println!("{:?}", features[i]);
        }
        
        let start = std::time::Instant::now();
        for i in 0..1 {
            let features = get_geometric_hash_as_u32_from_structure(&compact, HashType::PDBTrRosetta, 0, 0);
        }
        let duration = start.elapsed();
        println!("Non-SIMD feature computation: {:?}", duration);
        let features = get_geometric_hash_as_u32_from_structure(&compact, HashType::PDBTrRosetta, 0, 0);
        println!("Number of features: {}", features.len());
        for i in 0..10 {
            println!("{:?}", features[i]);
        }

    }
}