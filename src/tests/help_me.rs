#[derive(Default, Clone)]
struct Ends {
    n: i32,
}

trait EndExt {
    fn set_n(&mut self, n: i32);

    fn get_n(&self) -> i32;
}

impl EndExt for Ends {
    fn set_n(&mut self, n: i32) {
        self.n = n;
    }

    fn get_n(&self) -> i32 {
        self.n
    }
}

fn main() {
    let mut ends_vec = vec![Ends::default(); 100];

    do_something(ends_vec.as_mut_slice(), &Ends { n: 256 });
}

fn do_something<T: EndExt>(a: &mut [T], b: &T) {
    let mut sub_a = Vec::new();

    for e in a.iter_mut() {
        if e.get_n() % 2 == 0 {
            sub_a.push(e);
        }
    }

    // worker_in_deeper(sub_a.as_mut_slice(), b);
}

fn worker_in_deeper<T: EndExt>(a: &mut [T], b: &T) {
    // let mut s = String::new();
    // let sr = &mut s;

    // s.push_str("AS");
    // sr.push_str("AS");
    // s.push_str("AS");
    // sr.push_str("AS");
}
