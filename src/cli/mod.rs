// pub mod somemodule;
pub mod workflows;

#[derive(Debug)]
enum Subcommand {
    Index,
    Query,
    // Add subcommands here
}

pub enum AppArgs {
    Global {
        help: bool,
    },
    Index {
        pdb_dir: Option<String>,
        pdb_path_vec: Vec<String>,
        index_path: String,
        num_threads: usize,
        verbose: bool,
        help: bool,
    },
    Query {
        threads: usize,
        index_path: Option<String>,
        help: bool,
    },
}
