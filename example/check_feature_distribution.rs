use std::collections::HashMap;
// read tsv
use std::fs::File;
use std::io::{BufRead, BufReader, Write, Read};
use std::path::PathBuf;

use motifsearch::controller::{self, Controller, GeometryHashCollector};
use motifsearch::index::builder::IndexBuilder;
use motifsearch::index::IndexTablePrinter;
use motifsearch::PDBReader;

fn load_path(dir: &str) -> Vec<String> {
    // Load all pdbs in given path
    let mut pdb_paths = Vec::new();
    let paths = std::fs::read_dir(dir).expect("Unable to read pdb directory");
    for path in paths {
        let path = path.expect("Unable to read path");
        let path = path.path();
        let path = path.to_str().expect("Unable to convert path to string");
        // If the path is a pdb file, add it to the list
        if path.ends_with(".pdb") {
            pdb_paths.push(path.to_string());
        }
    }
    pdb_paths
}

#[derive(Clone, Debug)]
pub struct PeptidaseInfo {
    pub pdb_id: String,
    pub uniprot_id: String,
    pub subfamily: String,
    pub family: String,
    pub clan: String,
    pub name: String,
    pub active_site_raw: String,
    pub active_site: Vec<u64>,
    pub len_active_site: usize,
}
// Format for PeptidaseInfo
impl PeptidaseInfo {
    fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.pdb_id,
            self.uniprot_id,
            self.subfamily,
            self.family,
            self.clan,
            self.name,
            self.active_site_raw,
            self.active_site.len(),
        )
    }
}


fn load_info_tsv(path: &str) -> HashMap<String, PeptidaseInfo> {
    let mut info_tsv_reader = BufReader::new(File::open(path).expect("Failed opening file"));
    let mut info_tsv_lines = String::new();
    
    info_tsv_reader.read_to_string(&mut info_tsv_lines).expect("Failed reading line");
    let mut info_tsv_lines = info_tsv_lines.lines();
    let mut merops_info: HashMap<String, PeptidaseInfo> = HashMap::new();
    info_tsv_lines.next();
    for line in info_tsv_lines {
        let mut line = line.split("\t");
        let pdb_id = line.next().expect("Failed parsing line").to_string();
        let uniprot_id = line.next().expect("Failed parsing line").to_string();
        let subfamily = line.next().expect("Failed parsing line").to_string();
        let family = line.next().expect("Failed parsing line").to_string();
        let clan = line.next().expect("Failed parsing line").to_string();
        let name = line.next().expect("Failed parsing line").to_string();
        let active_site_raw = line.next().expect("Failed parsing line").to_string();
        // active site is a comma separated list of numbers starting with alphabet
        // Example: E275, E321, E418
        // ignore alphabet
        let active_site = active_site_raw
            .split(", ")
            .map(|x| x[1..]
            .parse::<u64>()
            .expect("Failed converting active site"))
            .collect::<Vec<u64>>();
        let len_active_site = active_site.len();
        let pdb_id_new = pdb_id.clone();
        let peptidase_info = PeptidaseInfo {
            pdb_id,
            uniprot_id,
            subfamily,
            family,
            clan,
            name,
            active_site_raw,
            active_site,
            len_active_site,
        };
        merops_info.insert(pdb_id_new, peptidase_info);
    }
    merops_info
}

fn calc_feature_for_only_active_site(pdb_path: &String, info: &PeptidaseInfo) ->
    Vec<(String, PeptidaseInfo, String, String, String, String, Vec<f32>, Vec<f32>)> 
{
    let pdb_reader = PDBReader::from_file(pdb_path).expect("pdb file not found");
    let structure = pdb_reader.read_structure().expect("structure read failed");
    let compact = structure.to_compact();
    let mut active_site_residue_serial = vec![];
    let mut active_site_residue_name = vec![];
    // iterate over residue_serial
    for (i, residue_serial) in compact.residue_serial.iter().enumerate() {
        if info.active_site.contains(residue_serial) {
            active_site_residue_serial.push(i);
            let residue_name = compact.residue_name[i].clone();
            active_site_residue_name.push(residue_name);
        }
    }
    // output: pdb_path, info, res1, res2, res1_name, res2_name, trr, pdb
    let mut output = vec![];
    for i in active_site_residue_serial.iter() {
        for j in active_site_residue_serial.iter() {
            if i == j {
                continue;
            }
            let trr = compact.get_trrosetta_feature(*i, *j).unwrap_or([0.0; 6]).to_vec();
            if trr[0] > 20.0 {
                continue;
            }
            let pdb_feature = compact.get_pdb_feature(*i, *j);
            let res1 = String::from_utf8(compact.residue_name[*i].clone().to_vec()).expect("Failed to convert residue name");
            let res2 = String::from_utf8(compact.residue_name[*j].clone().to_vec()).expect("Failed to convert residue name");
            let n =  compact.chain_per_residue[*i].to_string() + compact.residue_serial[*i].to_string().as_str();
            let m = compact.chain_per_residue[*j].to_string() + compact.residue_serial[*j].to_string().as_str();
            output.push((pdb_path.to_string(), info.clone(), n, m, res1, res2, trr, pdb_feature));
        }
    }
    output
}

fn main() {
    // Measure runtime
    let start = std::time::Instant::now();
    let pdb_paths = load_path("/Users/hbk/Projects/Lab/06_FoldMotif/repos/merops/merops_pdb/");
    let mut controller = Controller::new(pdb_paths);
    controller.fill_numeric_id_vec();
    controller.collect_hash();
    let end = std::time::Instant::now();
    println!("Elapsed time: {:?} for {:?} pdbs", end - start, controller.path_vec.len());
    // controller.save_raw_feature("data/merops_pdb_hash_raw.tsv", false);
    // controller.save_id_vec("data/merops_id_vec.tsv");

    // Read merops information
    let info_tsv_path = "/Users/hbk/Projects/Lab/06_FoldMotif/repos/merops/pdb_pa.tsv";
    let merops_info = load_info_tsv(info_tsv_path);
    
    // Calculate feature for only active site
    println!("Calculating feature for only active site");
    let mut output = vec![];
    for pdb_path in controller.path_vec.iter() {
        let pdb_path_str = pdb_path;
        let pdb_id = pdb_path.split("/").last().expect("Failed parsing name").split(".").next().expect("Failed parsing nam2");
        // println!("Processing {}", pdb_id);
        let info = merops_info.get(pdb_id);
        if info.is_none() {
            continue;
        }
        let output_for_pdb = calc_feature_for_only_active_site(pdb_path_str, info.unwrap());
        output.extend(output_for_pdb);
    }
    // Write output to tsv
    let mut output_tsv = File::create("data/merops_pdb_hash_active_site.tsv").expect("Failed creating output file");
    writeln!(output_tsv, "pdb_path\tpdb_id\tuniprot_id\tsubfamily\tfamily\tclan\tname\tactive_site_raw\tlen_active_site\t\tres1\tres2\tres1_name\tres2_name\ttrr_dist\ttrr_omega\ttrr_phi1\ttrr_phi2\ttrr_psi1\ttrr_psi2\tpdb_ca_dist\tpdb_cb_dist\tpdb_ca_cb_angle").expect("Failed writing header");
    for (pdb_path, info, res1, res2, res1_name, res2_name, trr, pdb_feature) in output {
        let output_line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            pdb_path, info.to_tsv(), res1, res2, res1_name, res2_name,
            trr[0], trr[1], trr[2], trr[3], trr[4], trr[5], pdb_feature[0], pdb_feature[1], pdb_feature[2]
        );
        writeln!(output_tsv, "{}", output_line).expect("Failed writing output");
    }
    
    // End
    println!("Done");
}
