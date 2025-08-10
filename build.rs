fn main() {
    if !cfg!(feature = "bench") {
        return;
    }
    println!("Running build.rs for benchmarks only...");
    cxx_build::bridge("src/cxx.rs")
        .file("source/lets_be_rational.cpp")
        .file("source/normaldistribution.cpp")
        .file("source/rationalcubic.cpp")
        .file("source/erf_cody.cpp")
        .flag("-finput-charset=UTF-8")
        .flag("-fextended-identifiers")
        .flag("-O3")
        .flag("-DNDEBUG")
        .flag("-ffp-contract=fast")
        .flag("-march=native")
        .flag("-Ofast")
        .flag("-flto")
        .flag("-w")
        .include("source")
        .flag_if_supported("-std=c++23")
        .compile("lets_be_rational");

    println!("cargo:rerun-if-changed=src/cxx.rs");
    println!("cargo:rerun-if-changed=source/lets_be_rational.h");
    println!("cargo:rerun-if-changed=source/lets_be_rational.cpp");
    println!("cargo:rerun-if-changed=source/erf_cody.cpp");
    println!("cargo:rerun-if-changed=source/normaldistribution.h");
    println!("cargo:rerun-if-changed=source/normaldistribution.cpp");
    println!("cargo:rerun-if-changed=source/rationalcubic.h");
    println!("cargo:rerun-if-changed=source/rationalcubic.cpp");
}
