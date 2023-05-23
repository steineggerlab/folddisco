use motifsearch::*;

fn main() {
    // Load directory
    let yeast_pdb_paths = motifsearch::test::load_path("data/yeast_full");
    let mut controller = motifsearch::controller::Controller::new(yeast_pdb_paths);
    controller.fill_numeric_id_vec();
    controller.save_raw_feature("yeast_raw_feature.tsv");
}