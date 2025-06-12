#[cfg(feature = "foldcomp")] 
extern crate cmake;

fn main() {
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
