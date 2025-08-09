fn main() {
    cxx_build::bridge("src/cxx.rs")
        .file("source/lets_be_rational.cpp")
        .file("source/normaldistribution.cpp")
        .file("source/rationalcubic.cpp")
        .file("source/erf_cody.cpp")
        .flag("-O3")
        .flag("-DNDEBUG")
        .flag("-fextended-identifiers")
        .flag("-ffp-contract=fast")
        .flag("-march=native")
        .include("source")
        .flag_if_supported("-std=c++17")
        .compile("lets_be_rational");

    println!("cargo:rerun-if-changed=src/cxx.rs");
    println!("cargo:rerun-if-changed=source/lets_be_rational.h");
    println!("cargo:rerun-if-changed=source/erf_cody.cpp");
    println!("cargo:rerun-if-changed=source/normaldistribution.h");
    println!("cargo:rerun-if-changed=source/rationalcubic.h");
}
