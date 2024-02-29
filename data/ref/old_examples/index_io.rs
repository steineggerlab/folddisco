
use std::{sync::{atomic::{AtomicUsize, Ordering}, Arc, Mutex}, time::Instant, io::{BufWriter, Write}, fs::File};

use byteorder::LittleEndian;
// This example is to check if the IndexAllocator works as intended
use motifsearch::{prelude::*, structure::atom::Atom};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

fn main() {
    const NUM_THREADS: usize = 6;
    
    // Load directory
    let yeast_pdb_paths = motifsearch::utils::loader::load_path("data/serine_peptidases_filtered");
    let mut controller = motifsearch::controller::Controller::new(yeast_pdb_paths);
    controller.num_threads_file = NUM_THREADS;
    controller.num_threads_hash = NUM_THREADS;
    controller.fill_numeric_id_vec();

    let mut id_lookup_map = FxHashMap::<String, usize>::default();
    for i in 0..controller.path_vec.len() {
        id_lookup_map.insert(controller.path_vec[i].clone(), controller.numeric_id_vec[i]);
    }
    let id_lookup_map = Arc::new(id_lookup_map);
    
    // Save offset
    let mut offset = AtomicUsize::new(0);
    let offset_vec = Arc::new(Mutex::new(
        Vec::<(u64, (usize, usize))>::with_capacity(controller.path_vec.len())
    ));
    
    rayon::ThreadPoolBuilder::new().num_threads(NUM_THREADS).build_global().unwrap();
    let start = Instant::now();
     let mut result = (&controller.path_vec).into_par_iter().map(
        |pdb_path| {
            // Read structure from PDB file
            let pdb_reader = PDBReader::from_file(&pdb_path).expect(
                &log_msg(FAIL, "Failed to read PDB file")
            );
            let compact = pdb_reader.read_structure().expect(
                &log_msg(FAIL, "Failed to read PDB file")
            );
            let compact = compact.to_compact();
            // Iterate over all residue pairs
            let res_bound = get_all_combination(compact.num_residues, false);
            let mut hash_collection: Vec<usize> = res_bound.into_par_iter().map(
               |(res1, res2)| {
                    if res1 == res2 {
                        return 0;
                    }
                    let feature = compact.get_trrosetta_feature(res1, res2).unwrap_or(
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                    );
                    if feature[0] < 2.0 || feature[0] > 20.0 {
                        return 0;
                    }
                    let hash = HashValue::perfect_hash(
                        feature[0], feature[1], feature[2], feature[3], feature[4], feature[5]
                    );
                    hash.as_usize()
                }
            ).collect();
            
            // Filter out invalid residue pairs & deduplicate
            // Remove zero
            hash_collection.sort_unstable();
            hash_collection.dedup();
            hash_collection.retain(|&x| x != 0);

            // Drop everything
            drop(compact);
            drop(pdb_reader);
            hash_collection.shrink_to_fit();
            // Check offset map is empty
            
            let offset = offset.fetch_add(hash_collection.len(), Ordering::Relaxed);

            // Append offset: numeric id as key, offset as value
            let nid = id_lookup_map.get(pdb_path.as_str()).unwrap().clone();
            let mut offset_vec = offset_vec.lock().unwrap();
            offset_vec.push((nid as u64, (offset, hash_collection.len())));

            hash_collection
        }
    ).flatten().collect::<Vec<_>>();
    result.shrink_to_fit();
    let result = unsafe {
        std::slice::from_raw_parts(
            result.as_ptr() as *const u8,
            result.len() * std::mem::size_of::<usize>()
        )
    };
    // Print size of result
    println!("Result size: {}", result.len());

    // Save to file
    let mut file = File::create("analysis/serine_peptidases.bin").expect(
        &log_msg(FAIL, "Failed to create index file")
    );
    let mut buf = BufWriter::new(file);
    buf.write_all(result).expect(
        &log_msg(FAIL, "Failed to write index file")
    );
    buf.flush().expect(
        &log_msg(FAIL, "Failed to flush index file")
    );

    let mut file = File::create("analysis/serine_peptidases.offset").expect(
        &log_msg(FAIL, "Failed to create offset file")
    );
    let mut writer = BufWriter::new(file);
    for (key, value) in offset_vec.lock().unwrap().iter() {
        writer.write_all(&key.to_le_bytes());
        writer.write_all(&value.0.to_le_bytes());
        writer.write_all(&value.1.to_le_bytes());
    }

    let end = Instant::now();
    println!("Time elapsed with {} threads: {:?}", NUM_THREADS, end - start);
    // Clear everything
}
