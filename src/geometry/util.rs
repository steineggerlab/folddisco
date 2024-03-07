// Discretizers 
pub fn discretize_f32_value_into_u64(val: f32, min: f32, max: f32, num_bin: f32) -> u64 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    ((val - min) * (disc_f) + 0.5) as u64
}

pub fn discretize_f32_value_into_u32(val: f32, min: f32, max: f32, num_bin: f32) -> u32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    ((val - min) * (disc_f) + 0.5) as u32
}

pub fn continuize_u64_value_into_f32(val: u64, min: f32, max: f32, num_bin: f32) -> f32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    (val as f32) * (cont_f) + min
}

pub fn continuize_u32_value_into_f32(val: u32, min: f32, max: f32, num_bin: f32) -> f32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    (val as f32) * (cont_f) + min
}

pub fn normalize_f32_value(val: f32, min: f32, max: f32) -> f32 {
    (val - min) / (max - min)
}


pub fn map_aa_to_u8(aa: &[u8; 3]) -> u8 {
    // Applied to handle the case of non-standard amino acids
    // reference: gemmi/blob/master/src/resinfo.cpp (https://github.com/project-gemmi/gemmi)
    match aa {
        b"ALA" | b"ABA" | b"ORN" | b"DAL" | b"AIB" | b"ALC" | b"MDO" | b"MAA" | b"DAB" => 0, // ALA, A, total 9
        b"ARG" | b"DAR" | b"CIR" | b"AGM" => 1, // ARG, R, total 4
        b"ASN" | b"DSG" | b"MEN" | b"SNN" => 2, // ASN, N, total 4
        b"ASP" | b"0TD" | b"DAS" | b"IAS" | b"PHD" | b"BFD" | b"ASX" => 3, // ASP, D, total 7, ASX is included here
        b"CYS" | b"CSO" | b"CSD" | b"CME" | b"OCS" | b"CAS" | b"CSX" | b"CSS" | 
        b"YCM" | b"DCY" | b"SMC" | b"SCH" | b"SCY" | b"CAF" | b"SNC" | b"SEC" => 4, // CYS, C, total 16, SEC is included here
        b"GLN" | b"DGN" | b"CRQ" | b"MEQ" => 5, // GLN, Q, total 4
        b"GLU" | b"PCA" | b"DGL" | b"CGU" | b"FGA" | b"B3E" | b"GLX" => 6, // GLU, E, total 7, GLX is included here
        b"GLY" | b"CR2" | b"SAR" | b"GHP" | b"GL3" => 7, // GLY, G, total 5
        b"HIS" | b"HIC" | b"DHI" | b"NEP" | b"CR8" | b"MHS" => 8, // HIS, H, total 6
        b"ILE" | b"DIL" => 9, // ILE, I, total 2
        b"LEU" | b"DLE" | b"NLE" | b"MLE" | b"MK8"=> 10, // LEU, L, total 5
        b"LYS" | b"KCX" | b"LLP" | b"MLY" | b"M3L" | b"ALY" | b"MLZ" | b"DLY" | 
        b"KPI" | b"PYL" => 11, // LYS, K, total 10, PYL is included here
        b"MET" | b"MSE" | b"FME" | b"NRQ" | b"CXM" | b"SME" | b"MHO" | b"MED" => 12, // MET, M, total 8
        b"PHE" | b"DPN" | b"PHI" | b"MEA" | b"PHL" => 13, // PHE, F, total 5
        b"PRO" | b"HYP" | b"DPR" => 14, // PRO, P, total 3
        b"SER" | b"CSH" | b"SEP" | b"DSN" | b"SAC" | b"GYS" | b"DHA" | b"OAS" => 15, // SER, S, total 8
        b"THR" | b"TPO" | b"CRO" | b"DTH" | b"BMT" | b"CRF" => 16, // THR, T, total 6
        b"TRP" | b"DTR" | b"TRQ" | b"TOX" | b"0AF" => 17, // TRP, W, total 5
        b"TYR" | b"PTR" | b"TYS" | b"TPQ" | b"DTY" | b"OMY" => 18,
        b"VAL" | b"DVA" | b"MVA" | b"FVA" => 19,
        _ => 255,
    }
}
pub fn map_u8_to_aa(aa: u8) -> &'static str {
    match aa {
        0 => "ALA",
        1 => "ARG",
        2 => "ASN",
        3 => "ASP",
        4 => "CYS",
        5 => "GLN",
        6 => "GLU",
        7 => "GLY",
        8 => "HIS",
        9 => "ILE",
        10 => "LEU",
        11 => "LYS",
        12 => "MET",
        13 => "PHE",
        14 => "PRO",
        15 => "SER",
        16 => "THR",
        17 => "TRP",
        18 => "TYR",
        19 => "VAL",
        _ => "UNK",
    }
}

pub fn map_aa_pair_to_u32(aa1: &[u8; 3], aa2: &[u8; 3]) -> u32 {
    let output = (map_aa_to_u8(aa1) as u32) * 20 + map_aa_to_u8(aa2) as u32;
    assert!(output < 512);
    output
}

pub fn map_aa_u32_pair_to_u32(aa1: u32, aa2: u32) -> u32 {
    let output = aa1 * 20 + aa2;
    assert!(output < 512);
    output
}

pub fn map_u32_to_aa_u32_pair(pair: u32) -> (u32, u32) {
    let aa1 = (pair / 20) as u32;
    let aa2 = (pair % 20) as u32;
    (aa1, aa2)
}

pub fn map_u32_to_aa_pair(pair: u32) -> (String, String) {
    let aa1 = (pair / 20) as u8;
    let aa2 = (pair % 20) as u8;
    (map_u8_to_aa(aa1).to_string() , map_u8_to_aa(aa2).to_string())
}