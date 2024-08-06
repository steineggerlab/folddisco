
// #[cfg(feature = "foldcomp")] 
// extern crate bindgen;
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
    // println!("cargo:rustc-link-lib=c++"); // TODO: FIXME

    // Determine the platform and link the appropriate C++ standard library
    if cfg!(target_os = "windows") {
        println!("cargo:rustc-link-lib=dylib=foldcomp_ffi");
        println!("cargo:rustc-link-lib=dylib=msvcrt");
    } else if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=static=foldcomp_ffi");
        println!("cargo:rustc-link-lib=dylib=c++");
    } else if cfg!(target_os = "linux") {
        // Detect the C++ compiler used by CMake and link the appropriate standard library
        let cxx_compiler = std::env::var("CXX").unwrap_or_else(|_| "c++".to_string());
        if cxx_compiler.contains("clang") {
            println!("cargo:rustc-link-lib=static=foldcomp_ffi");
            println!("cargo:rustc-link-lib=dylib=c++");
        } else {
            println!("cargo:rustc-link-lib=static=foldcomp_ffi");
            println!("cargo:rustc-link-lib=dylib=stdc++");
        }
    }

    // Generate bindings. If regeneration is needed, uncomment the following code
    // let bindings = bindgen::Builder::default()
    //     .header("lib/foldcomp/foldcompffi.h")
    //     .generate()
    //     .expect("Unable to generate bindings");

    // let out_path = PathBuf::from("lib/foldcomp");
    // bindings
    //     .write_to_file(out_path.join("bindings.rs"))
    //     .expect("Couldn't write bindings!");

    }
}