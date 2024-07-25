
#[cfg(feature = "foldcomp")] 
extern crate bindgen;
#[cfg(feature = "foldcomp")] 
extern crate cmake;

use std::env;
use std::path::PathBuf;

fn main() {
    #[cfg(feature = "foldcomp")] 
    {
    let dst = cmake::Config::new("lib/foldcomp")
        .define("BUILD_FFI", "ON")
        .build_target("foldcomp_ffi")
        .build();

    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-lib=static=foldcomp_ffi");
    println!("cargo:rustc-link-lib=c++"); // TODO: FIXME

    let bindings = bindgen::Builder::default()
        .header("lib/foldcomp/foldcompffi.h")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from("lib/foldcomp");
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");

    }
}