// Time load_lookup_from_file, cold cache then warm.
use folddisco::index::lookup::load_lookup_from_file;

fn main() {
    let path = std::env::args().nth(1).expect("usage: lookup_load_bench <lookup>");
    let t = std::time::Instant::now();
    let table = load_lookup_from_file(&path);
    let elapsed = t.elapsed();
    // Touch a few entries so the timing cannot be hiding deferred work.
    let n = table.len();
    let probe: Vec<String> = [0usize, n / 3, n / 2, n - 1].iter()
        .filter(|&&i| i < n)
        .map(|&i| { let e = table.entry(i); format!("{} id={} nres={} plddt={} key={}", e.name, e.id, e.nres, e.plddt, e.db_key) })
        .collect();
    println!("entries      : {}", n);
    println!("load          : {:?}", elapsed);
    for p in probe { println!("  probe       : {}", p); }
    if std::env::var("SCAN_ALL").is_ok() {
        let t2 = std::time::Instant::now();
        let total: u64 = table.records().iter().map(|r| r.nres).sum();
        println!("full scan    : {:?} (nres sum {})", t2.elapsed(), total);
    }
}
