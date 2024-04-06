use proc_macro2::{Span, TokenStream};
use quote::quote;
use syn::{parse::Parse, parse_macro_input};

#[proc_macro]
/// # Example
///
/// ```
/// set_mlog!(stringify!(MarkDuplicates))
/// ```
pub fn set_mlog(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let minput = parse_macro_input!(input as MInput);

    let q = make_mlog_mod(&minput);

    eprintln!("{}", q);

    q.into()
}

const DEBUG: &'static str = "debug";
const INFO: &'static str = "info";
const WARN: &'static str = "warn";
const ERROR: &'static str = "error";
fn make_mlog_mod(minput: &MInput) -> TokenStream {
    let log_macro_ts = [DEBUG, INFO, WARN, ERROR]
        .into_iter()
        .map(|level| define_log_macro(minput, level));

    quote!(
        mod mlog {
            const EXIST:&'static str = "HEY!";

            #(
                #log_macro_ts
            )*
        }
    )
}

fn define_log_macro(minput: &MInput, level: &str) -> TokenStream {
    let level = syn::Ident::new(level, Span::call_site());
    let target = &minput.target;
   
    let macro_ident = if level == WARN {
        syn::Ident::new("warns", Span::call_site())
    } else {
        level.clone()
    };

    quote!(
        //#[macro_export]
        #[allow(unused_macros)]
        macro_rules! #macro_ident {
            ($($tt:tt)+) => {
                log::#level!(target: #target, $($tt)+)
            }
        }

        pub(super) use #macro_ident as #level;
    )
}

struct MInput {
    target: syn::Expr,
}

impl Parse for MInput {
    fn parse(input: syn::parse::ParseStream) -> syn::Result<Self> {
        Ok(Self {
            target: input.parse()?,
        })
    }
}
