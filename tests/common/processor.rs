use motifsearch::structure::feature::Torsion;
use motifsearch::structure::core;
use motifsearch::PDBReader;

#[derive(Debug)]
pub struct TmpGeometry {
    residue_pair: Vec<(u64, u64)>,
    angle: Vec<f32>,
    distance: Vec<f32>,

    torsion: Torsion,
}

impl TmpGeometry {
    fn new() -> TmpGeometry {
        TmpGeometry {
            residue_pair: Vec::new(),
            angle: Vec::new(),
            distance: Vec::new(),

            torsion: Torsion::new(),
        }
    }

    fn push(&mut self, res_pair: (u64, u64), angle: f32, distance: f32) {
        self.residue_pair.push(res_pair);
        self.angle.push(angle);
        self.distance.push(distance);
    }

    fn push_torsion(&mut self, torsion: Torsion) {
        self.torsion = torsion;
    }
}

pub struct TmpPDB {
    //IDEA: replace Vec to HashMap
    pdbs: Vec<String>,
    geometrys: Vec<TmpGeometry>,
}

impl TmpPDB {
    fn new() -> TmpPDB {
        TmpPDB {
            pdbs: Vec::new(),
            geometrys: Vec::new(),
        }
    }

    fn push(&mut self, pdb: String, geometry: TmpGeometry) {
        self.pdbs.push(pdb);
        self.geometrys.push(geometry);
    }
}

pub fn process_pdbs(path_vec: &Vec<String>) -> TmpPDB {
    /* to get distributions of dist, angle, torsion(psi) angle
    similar to controllers but not hash yet */
    let mut tmp_pdb = TmpPDB::new();

    for i in 0..path_vec.len() {
        let pdb_path = &path_vec[i];
        let pdb_reader = PDBReader::from_file(pdb_path).expect("pdb file not found");
        let structure = pdb_reader.read_structure().expect("structure read failed");
        let compact = structure.to_compact();

        let torsion = structure.get_torsion();

        let mut tmp_geometry = process_structure(&compact);
        tmp_geometry.push_torsion(torsion);

        tmp_pdb.push(pdb_path.to_string(), tmp_geometry);
    }
    tmp_pdb
}

pub fn process_structure(compact: &core::CompactStructure) -> TmpGeometry {
    let mut tmp_geometry = TmpGeometry::new();
    /* process one pdb and return TmpGeometry */
    for i in 0..compact.num_residues {
        for j in i + 1..compact.num_residues {
            let resi_pair = compact.get_res_serial(i, j);
            let dist = compact
                .get_distance(i, j)
                .expect("compact failed to get distance");
            let angle = compact.get_angle(i, j).unwrap_or(0.0);
            // let angle = compact.get_angle(i,j).expect("compact failed to get angle");
            tmp_geometry.push(resi_pair, angle, dist);
        }
    }
    tmp_geometry
}
