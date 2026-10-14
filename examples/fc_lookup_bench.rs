// Time FoldcompLookup::load, cold cache then warm.
#[cfg(feature = "foldcomp")]
fn main() {
    use folddisco::structure::io::fcz::FoldcompLookup;
    let db = std::env::args().nth(1).expect("usage: fc_lookup_bench <db prefix>");
    let t = std::time::Instant::now();
    let lookup = FoldcompLookup::load(&db).expect("load failed");
    let elapsed = t.elapsed();
    println!("entries      : {}", lookup.len());
    println!("load         : {:?}", elapsed);
    let n = lookup.len();
    for i in [0usize, n / 2, n - 1] {
        if i < n { println!("  probe      : key={} name={}", lookup.key(i), lookup.name(i)); }
    }
    // key -> name, the per-hit shape
    let t2 = std::time::Instant::now();
    let probe_key = lookup.key(n / 3);
    let name = lookup.name_of_key(probe_key);
    println!("name_of_key  : {:?} -> {:?}", t2.elapsed(), name);
    // name -> key, the once-per-query shape. Opt-in, because it scans every
    // name and so dominates the process's resident set.
    if std::env::var("NAME_SCAN").is_ok() {
        let wanted = lookup.name(n - 7).to_string();
        let t3 = std::time::Instant::now();
        let key = lookup.key_of_name(&wanted);
        println!("key_of_name  : {:?} -> {:?} for {}", t3.elapsed(), key, wanted);
    }
}

#[cfg(not(feature = "foldcomp"))]
fn main() { eprintln!("needs the foldcomp feature"); }
