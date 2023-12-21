
use std::{sync::{atomic::{AtomicUsize, Ordering}, Arc}, time::Instant, io::{BufWriter, Write}, fs::File};

use byteorder::LittleEndian;
// This example is to check if the IndexAllocator works as intended
use motifsearch::{prelude::*, structure::atom::Atom, index::alloc::resize_allocation};
use rayon::prelude::*;

fn main() {
    const NUM_THREADS: usize = 6;
    
    // Load directory
    let yeast_pdb_paths = motifsearch::utils::loader::load_path("data/serine_peptidases_filtered");
    let mut controller = motifsearch::controller::Controller::new(yeast_pdb_paths);
    controller.num_threads_file = NUM_THREADS;
    controller.num_threads_hash = NUM_THREADS;
    controller.fill_numeric_id_vec();
    
    let alloc_size = 
    
    let mut allocator = IndexAllocator::new(controller.num_threads_hash, alloc_size);
    let mut dedup_size: AtomicUsize = AtomicUsize::new(0);
    
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
            // Single
            let mut hash_collection: Vec<usize> = res_bound.into_par_iter().enumerate().map(
                |(i, (res1, res2))| {
                        if res1 == res2 {
                            return 0usize;
                        }
                        let feature = compact.get_trrosetta_feature(res1, res2).unwrap_or(
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        );
                        if feature[0] < 2.0 || feature[0] > 20.0 {
                            return 0usize;
                        }
                        let hash = HashValue::perfect_hash(
                            feature[0], feature[1], feature[2], feature[3], feature[4], feature[5]
                        );
                        hash.as_usize()
                    }
            ).collect::<Vec<_>>();
            // Filter out invalid residue pairs & deduplicate
            // Remove zero
            hash_collection.sort_unstable();
            hash_collection.dedup();
            hash_collection.retain(|&x| x != 0);
            hash_collection.shrink_to_fit();
            
            // Check if reallocation is needed
            let len = hash_collection.len();
            let cursor = allocator.cursor.load(Ordering::Relaxed);
            let size = allocator.data_size.load(Ordering::Relaxed);
            if cursor + len >= size {
                // Use unsafe code to resize the allocation vector
                // Allow mutable method to be called in here
                unsafe {
                    allocator.data_size.store(size * 2, Ordering::Relaxed);
                    let size_as_mb = size as f32 * 8.0 / 1024.0 / 1024.0;
                    println!("Resized allocation vector to {}MB", size * 2);
                    // Get mutable reference to allocator with raw pointer
                    let mut ptr = allocator.allocation.as_ptr() as &mut [AtomicUsize];
                    // ptr.set_len(size * 2);
                    // Resize the allocation vector
                    resize_allocation(ptr, size * 2);
                    // Update the allocation vector
                }
            }
            
            // Allocate
            // allocator.fill_and_update(hash_collection);
            allocator.fill_usize_vec(hash_collection);
            
            // allocator.fill_usize_vec(hash_collection);
            
        }
    ).collect::<Vec<_>>();

    /** 
     * IMPLEMENTATION WITH VECTOR START
     *     let mut hash_collection: Vec<usize> = res_bound.into_par_iter().map(
     *           |(res1, res2)| {
     *                if res1 == res2 {
     *                    return 0;
     *                }
     *                let feature = compact.get_trrosetta_feature(res1, res2).unwrap_or(
     *                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
     *                );
     *                if feature[0] < 2.0 || feature[0] > 20.0 {
     *                    return 0;
     *                }
     *                let hash = HashValue::perfect_hash(
     *                    feature[0], feature[1], feature[2], feature[3], feature[4], feature[5]
     *                );
     *                hash.as_usize()
     *            }
     *        ).collect();
     *        
     *        // Filter out invalid residue pairs & deduplicate
     *        // Remove zero
     *        hash_collection.sort_unstable();
     *        hash_collection.dedup();
     *        hash_collection.retain(|&x| x != 0);
     *        // allocator.fill_and_update(hash_collection);
     *        // allocator.fill_usize_vec(hash_collection);
     *        // Drop everything
     *        drop(compact);
     *        drop(pdb_reader);
     *        hash_collection.shrink_to_fit();
     *        hash_collection
     *    }
     *).flatten().collect::<Vec<_>>();
     *result.shrink_to_fit();
     *let result = unsafe {
     *    std::slice::from_raw_parts(
     *        result.as_ptr() as *const u8,
     *        result.len() * std::mem::size_of::<usize>()
     *    )
     *};
     *   // Print size of result
     *   println!("Result size: {}", result.len());
     */

    // Resize the allocation
    allocator.save_to_file("analysis/raw_ecoli.bin").expect(
        &log_msg(FAIL, "Failed to save index file")
    );
    // // Save to file
    // let mut file = File::create("analysis/raw_ecoli.bin").expect(
    //     &log_msg(FAIL, "Failed to create index file")
    // );
    // let mut buf = BufWriter::new(file);
    // buf.write_all(result).expect(
    //     &log_msg(FAIL, "Failed to write index file")
    // );
    // buf.flush().expect(
    //     &log_msg(FAIL, "Failed to flush index file")
    // );
    
    let end = Instant::now();
    println!("Time elapsed with {} threads: {:?}", NUM_THREADS, end - start);
    // Clear everything
}
