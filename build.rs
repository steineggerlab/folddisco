#[cfg(feature = "foldcomp")] 
extern crate cmake;

use std::process::Command;

fn main() {
    // Set git version 
    // First read environment variable - FOLDDISCO_BUILD_VERSION
    let git_hash = std::env::var("FOLDDISCO_BUILD_VERSION").unwrap_or_else(|_| {
        let git_hash = Command::new("git").args(&["rev-parse", "HEAD"])
            .output()
            .ok()
            .and_then(|output| {
                if output.status.success() {
                    Some(String::from_utf8_lossy(&output.stdout).trim().to_string())
                } else {
                    None
                }
            });
        git_hash.unwrap_or_else(|| "unknown".into())
    });
    
    // Set environment variable
    println!("cargo:rustc-env=FOLDDISCO_BUILD_VERSION={}", git_hash);
    // Configure rebuild triggers
    println!("cargo:rerun-if-env-changed=FOLDDISCO_BUILD_VERSION");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=.git/HEAD");        

    // Build Foldcomp C++ library needed for Foldcomp support
    #[cfg(feature = "foldcomp")] 
    {
    let dst = cmake::Config::new("lib/foldcomp")
        .define("BUILD_FFI", "ON")
        .build_target("foldcomp_ffi")
        .build();

    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-lib=static=foldcomp_ffi");

    // Determine the platform and link the appropriate C++ standard library
    if cfg!(target_os = "windows") {
        println!("cargo:rustc-link-lib=dylib=foldcomp_ffi");
        println!("cargo:rustc-link-lib=dylib=msvcrt");
    } else if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=static=foldcomp_ffi");
        println!("cargo:rustc-link-lib=dylib=c++");
    } else if cfg!(target_os = "linux") {
        println!("cargo:rustc-link-lib=static=foldcomp_ffi");
        let link_type = std::env::var("FOLDCOMP_STD_LINK")
            .unwrap_or_else(|_| "dylib".into());
        let stdlib = std::env::var("FOLDCOMP_STD_LIB").unwrap_or_else(|_| {
            // Detect the C++ compiler used by CMake and link the appropriate standard library
            let cxx = std::env::var("CXX").unwrap_or_else(|_| "c++".into());
            if cxx.contains("clang") {
                "c++".into()
            } else {
                "stdc++".into()
            }
        });
        println!("cargo:rustc-link-lib={}={}", link_type, stdlib);
        if let Ok(dir) = std::env::var("FOLDCOMP_LIB_DIR") {
            println!("cargo:rustc-link-search=native={}", dir);
        }
    }

    }
}
