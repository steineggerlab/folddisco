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
    match aa {
        b"ALA" => 0,
        b"ARG" => 1,
        b"ASN" => 2,
        b"ASP" => 3,
        b"CYS" => 4,
        b"GLN" => 5,
        b"GLU" => 6,
        b"GLY" => 7,
        b"HIS" => 8,
        b"ILE" => 9,
        b"LEU" => 10,
        b"LYS" => 11,
        b"MET" => 12,
        b"PHE" => 13,
        b"PRO" => 14,
        b"SER" => 15,
        b"THR" => 16,
        b"TRP" => 17,
        b"TYR" => 18,
        b"VAL" => 19,
        _ => panic!("Invalid AA"),
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
        _ => panic!("Invalid AA"),
    }
}
