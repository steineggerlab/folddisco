
use std::collections::HashMap;
use lazy_static::lazy_static;

// Constants
// 1. for cb_dist
pub const MIN_DIST: f32 = 2.0;
pub const MAX_DIST: f32 = 20.0;
pub const NBIN_DIST: f32 = 8.0;
// 2. NEW IDEA for encoding angles; represent as sin and cos
pub const MIN_SIN_COS: f32 = -1.0;
pub const MAX_SIN_COS: f32 = 1.0;
pub const NBIN_TORSION_SIN_COS: f32 = 3.0;
pub const NBIN_PLANE_SIN_COS: f32 = 3.0;
pub const NBIN_SIN_COS: f32 = 3.0;
// Bitmasks
pub const BITMASK32_2BIT: u32 = 0x00000003;
pub const BITMASK32_3BIT: u32 = 0x00000007;
pub const BITMASK32_4BIT: u32 = 0x0000000F;
pub const BITMASK32_5BIT: u32 = 0x0000001F;
pub const BITMASK32_9BIT: u32 = 0x000001FF;

pub const BITMASK64_4BIT: u64 = 0x000000000000000F;
pub const BITMASK64_5BIT: u64 = 0x000000000000001F;

// Discretizers 
#[inline(always)]
pub fn discretize_f32_value_into_u64(val: f32, min: f32, max: f32, num_bin: f32) -> u64 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    ((val - min) * (disc_f) + 0.5) as u64
}

#[inline(always)]
pub fn discretize_f32_value_into_u32(val: f32, min: f32, max: f32, num_bin: f32) -> u32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    ((val - min) * (disc_f) + 0.5) as u32
}
#[inline(always)]
pub fn continuize_u64_value_into_f32(val: u64, min: f32, max: f32, num_bin: f32) -> f32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    (val as f32) * (cont_f) + min
}
#[inline(always)]
pub fn continuize_u32_value_into_f32(val: u32, min: f32, max: f32, num_bin: f32) -> f32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    (val as f32) * (cont_f) + min
}
// #[inline(always)]
pub fn normalize_f32_value(val: f32, min: f32, max: f32) -> f32 {
    (val - min) / (max - min)
}


// pub fn map_aa_to_u8(aa: &[u8; 3]) -> u8 {
//     // Applied to handle the case of non-standard amino acids
//     // reference: gemmi/blob/master/src/resinfo.cpp (https://github.com/project-gemmi/gemmi)
//     match aa {
//         b"ALA" | b"ABA" | b"ORN" | b"DAL" | b"AIB" | b"ALC" | b"MDO" | b"MAA" | b"DAB" => 0, // ALA, A, total 9
//         b"ARG" | b"DAR" | b"CIR" | b"AGM" => 1, // ARG, R, total 4
//         b"ASN" | b"DSG" | b"MEN" | b"SNN" => 2, // ASN, N, total 4
//         b"ASP" | b"0TD" | b"DAS" | b"IAS" | b"PHD" | b"BFD" | b"ASX" => 3, // ASP, D, total 7, ASX is included here
//         b"CYS" | b"CSO" | b"CSD" | b"CME" | b"OCS" | b"CAS" | b"CSX" | b"CSS" | 
//         b"YCM" | b"DCY" | b"SMC" | b"SCH" | b"SCY" | b"CAF" | b"SNC" | b"SEC" => 4, // CYS, C, total 16, SEC is included here
//         b"GLN" | b"DGN" | b"CRQ" | b"MEQ" => 5, // GLN, Q, total 4
//         b"GLU" | b"PCA" | b"DGL" | b"CGU" | b"FGA" | b"B3E" | b"GLX" => 6, // GLU, E, total 7, GLX is included here
//         b"GLY" | b"CR2" | b"SAR" | b"GHP" | b"GL3" => 7, // GLY, G, total 5
//         b"HIS" | b"HIC" | b"DHI" | b"NEP" | b"CR8" | b"MHS" => 8, // HIS, H, total 6
//         b"ILE" | b"DIL" => 9, // ILE, I, total 2
//         b"LEU" | b"DLE" | b"NLE" | b"MLE" | b"MK8"=> 10, // LEU, L, total 5
//         b"LYS" | b"KCX" | b"LLP" | b"MLY" | b"M3L" | b"ALY" | b"MLZ" | b"DLY" | 
//         b"KPI" | b"PYL" => 11, // LYS, K, total 10, PYL is included here
//         b"MET" | b"MSE" | b"FME" | b"NRQ" | b"CXM" | b"SME" | b"MHO" | b"MED" => 12, // MET, M, total 8
//         b"PHE" | b"DPN" | b"PHI" | b"MEA" | b"PHL" => 13, // PHE, F, total 5
//         b"PRO" | b"HYP" | b"DPR" => 14, // PRO, P, total 3
//         b"SER" | b"CSH" | b"SEP" | b"DSN" | b"SAC" | b"GYS" | b"DHA" | b"OAS" => 15, // SER, S, total 8
//         b"THR" | b"TPO" | b"CRO" | b"DTH" | b"BMT" | b"CRF" => 16, // THR, T, total 6
//         b"TRP" | b"DTR" | b"TRQ" | b"TOX" | b"0AF" => 17, // TRP, W, total 5
//         b"TYR" | b"PTR" | b"TYS" | b"TPQ" | b"DTY" | b"OMY" => 18,
//         b"VAL" | b"DVA" | b"MVA" | b"FVA" => 19,
//         _ => 255,
//     }
// }

lazy_static! {
    static ref AA_MAP: HashMap<&'static [u8; 3], u8> = {
        let mut m: HashMap<&'static [u8; 3], u8> = HashMap::new();
        // ALA group
        let ala_codes = [b"ALA", b"ABA", b"ORN", b"DAL", b"AIB", b"ALC", b"MDO", b"MAA", b"DAB"];
        for code in ala_codes.iter() {
            m.insert(code, 0);
        }
        // ARG group
        let arg_codes = [b"ARG", b"DAR", b"CIR", b"AGM"];
        for code in arg_codes.iter() { m.insert(code, 1); }
        // ASN group
        let asn_codes = [b"ASN", b"DSG", b"MEN", b"SNN"];
        for code in asn_codes.iter() { m.insert(code, 2); }
        // ASP group
        let asp_codes = [b"ASP", b"0TD", b"DAS", b"IAS", b"PHD", b"BFD", b"ASX"];
        for code in asp_codes.iter() { m.insert(code, 3); }
        // CYS group (total 16 including SEC)
        let cys_codes = [b"CYS", b"CSO", b"CSD", b"CME", b"OCS", b"CAS", b"CSX", b"CSS", b"YCM", b"DCY", b"SMC", b"SCH", b"SCY", b"CAF", b"SNC", b"SEC"];
        for code in cys_codes.iter() {
            m.insert(code, 4);
        }
        // GLN group
        let gln_codes = [b"GLN", b"DGN", b"CRQ", b"MEQ"];
        for code in gln_codes.iter() { m.insert(code, 5); }
        // GLU group
        let glu_codes = [b"GLU", b"PCA", b"DGL", b"CGU", b"FGA", b"B3E", b"GLX"];
        for code in glu_codes.iter() { m.insert(code, 6); }
        // GLY group
        let gly_codes = [b"GLY", b"CR2", b"SAR", b"GHP", b"GL3"];
        for code in gly_codes.iter() { m.insert(code, 7); }
        // HIS group
        let his_codes = [b"HIS", b"HIC", b"DHI", b"NEP", b"CR8", b"MHS"];
        for code in his_codes.iter() { m.insert(code, 8); }
        // ILE group
        let ile_codes = [b"ILE", b"DIL"];
        for code in ile_codes.iter() { m.insert(code, 9); }
        // LEU group
        let leu_codes = [b"LEU", b"DLE", b"NLE", b"MLE", b"MK8"];
        for code in leu_codes.iter() { m.insert(code, 10); }
        // LYS group
        let lys_codes = [b"LYS", b"KCX", b"LLP", b"MLY", b"M3L", b"ALY", b"MLZ", b"DLY", b"KPI", b"PYL"];
        for code in lys_codes.iter() { m.insert(code, 11); }
        // MET group
        let met_codes = [b"MET", b"MSE", b"FME", b"NRQ", b"CXM", b"SME", b"MHO", b"MED"];
        for code in met_codes.iter() { m.insert(code, 12); }
        // PHE group
        let phe_codes = [b"PHE", b"DPN", b"PHI", b"MEA", b"PHL"];
        for code in phe_codes.iter() { m.insert(code, 13); }
        // PRO group
        let pro_codes = [b"PRO", b"HYP", b"DPR"];
        for code in pro_codes.iter() { m.insert(code, 14); }
        // SER group
        let ser_codes = [b"SER", b"CSH", b"SEP", b"DSN", b"SAC", b"GYS", b"DHA", b"OAS"];
        for code in ser_codes.iter() { m.insert(code, 15); }
        // THR group
        let thr_codes = [b"THR", b"TPO", b"CRO", b"DTH", b"BMT", b"CRF"];
        for code in thr_codes.iter() { m.insert(code, 16); }
        // TRP group
        let trp_codes = [b"TRP", b"DTR", b"TRQ", b"TOX", b"0AF"];
        for code in trp_codes.iter() { m.insert(code, 17); }
        // TYR group
        let tyr_codes = [b"TYR", b"PTR", b"TYS", b"TPQ", b"DTY", b"OMY"];
        for code in tyr_codes.iter() { m.insert(code, 18); }
        // VAL group
        let val_codes = [b"VAL", b"DVA", b"MVA", b"FVA"];
        for code in val_codes.iter() { m.insert(code, 19); }
        m
    };
}

pub fn map_aa_to_u8(aa: &[u8; 3]) -> u8 {
    *AA_MAP.get(aa).unwrap_or(&255)
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

#[inline(always)]
pub fn map_aa_pair_to_u32(aa1: &[u8; 3], aa2: &[u8; 3]) -> u32 {
    let output = (map_aa_to_u8(aa1) as u32) * 20 + map_aa_to_u8(aa2) as u32;
    assert!(output < 512);
    output
}

#[inline(always)]
pub fn map_aa_u32_pair_to_u32(aa1: u32, aa2: u32) -> u32 {
    let output = aa1 * 20 + aa2;
    assert!(output < 512);
    output
}

#[inline(always)]
pub fn map_u32_to_aa_u32_pair(pair: u32) -> (u32, u32) {
    let aa1 = (pair / 20) as u32;
    let aa2 = (pair % 20) as u32;
    (aa1, aa2)
}

#[inline(always)]
pub fn map_u32_to_aa_pair(pair: u32) -> (String, String) {
    let aa1 = (pair / 20) as u8;
    let aa2 = (pair % 20) as u8;
    (map_u8_to_aa(aa1).to_string() , map_u8_to_aa(aa2).to_string())
}

pub fn map_one_letter_to_u8_vec(aa: char) -> Vec<u8> {
    match aa {
        'A' => vec![0],
        'R' => vec![1],
        'N' => vec![2],
        'D' => vec![3],
        'C' => vec![4],
        'Q' => vec![5],
        'E' => vec![6],
        'G' => vec![7],
        'H' => vec![8],
        'I' => vec![9],
        'L' => vec![10],
        'K' => vec![11],
        'M' => vec![12],
        'F' => vec![13],
        'P' => vec![14],
        'S' => vec![15],
        'T' => vec![16],
        'W' => vec![17],
        'Y' => vec![18],
        'V' => vec![19],
        // Handle the case of non-standard amino acids
        'B' => vec![2, 3], // Asn, Asp
        'Z' => vec![5, 6], // Gln, Glu
        'X' => vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], // All
        'J' => vec![9, 10], // Ile, Leu
        'U' => vec![4], // Selenocysteine 
        'O' => vec![11], // Pyrrolysine
        // Custom characters to represent groups of amino acids
        // Positively charged: +, Negatively charged: -, Polar or hydrophilic: ^, Non-polar or hydrophobic: ~, Aromatic: @
        '+' => vec![1, 8, 11], // Arg, His, Lys
        '-' => vec![3, 6], // Asp, Glu
        '^' => vec![2, 5, 15, 16, 18], // Asn, Gln, Ser, Thr, Tyr
        '~' => vec![0, 7, 9, 10, 12, 13, 14, 19], // Ala, Gly, Ile, Leu, Met, Phe, Pro, Val
        '@' => vec![8, 13, 17, 18], // His, Phe, Trp, Tyr
        _ => vec![255], // Unknown
    }
}

#[inline(always)]
pub fn is_aa_group_char(c: char) -> bool {
    if c.is_ascii_alphabetic() {
        true
    } else {
        match c {
            '+' | '-' | '^' | '~' | '@' => true,
             _ => false,
        }
    }
}
