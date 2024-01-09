mod normal_distribution;
mod erf_cody;

#[test]
fn main() {
    let x = erf_cody::calerf(0.01, 1);
    println!("{x}");
    let x = erf_cody::calerf(1.0, 1);
    println!("{x}");
    let x = erf_cody::calerf(4.0, 1);
    println!("{x}");
    let x = erf_cody::calerf(4.000000001, 1);
    println!("{x}");
    let x = erf_cody::calerf(4.0000001, 1);
    println!("{x}");
    let x = erf_cody::calerf(4.00001, 1);
    println!("{x}");
}