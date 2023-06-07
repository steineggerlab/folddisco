use crate::cli::*;

// TODO: Generated. need to be

const HELP_INDEX: &str = "\
USAGE: motifsearch index [OPTIONS] <PDBS...>
Options:
    -d, --pdb-dir <PDB_DIR>     Directory containing PDB files
    -i, --index-path <INDEX_PATH>   Path to save the index table
    -t, --threads <THREADS>     Number of threads to use
    -v, --verbose               Print verbose messages
    -h, --help                  Print this help menu
";

pub fn build_index(env: AppArgs) {
    match env {
        AppArgs::Index {
            pdb_dir,
            pdb_path_vec,
            index_path,
            num_threads,
            verbose,
            help,
        } => {
            if help {
                println!("{}", HELP_INDEX);
            } else {
                if verbose {
                    println!("Building index table...");
                }
                todo!("IMPLEMENT INDEX BUILDING");
                // let mut index = Index::new();
                // if let Some(pdb_dir) = pdb_dir {
                //     index.build_from_dir(&pdb_dir, num_threads, verbose);
                // } else {
                //     index.build_from_files(&pdb_path_vec, num_threads, verbose);
                // }
                // index.save(&index_path);
                if verbose {
                    println!("Done!");
                }
            }
        }
        _ => {
            println!("Invalid subcommand. Try `motifsearch --help` for more information.");
        }
    }
}
