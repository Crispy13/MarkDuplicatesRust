use macro_sup::set_mlog;


set_mlog!("TOAST!");

fn log_test() {
    mlog::error!("HEY!");
}

fn main() {
    mlog::warn!("HEY!");
}