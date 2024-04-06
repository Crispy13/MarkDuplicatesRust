use macro_sup::set_mlog;

const s: &str = stringify!(TOAST);
set_mlog!(s);

fn log_test() {
    mlog::error!("HEY!");
}

fn main() {
    mlog::warn!("HEY!");
}