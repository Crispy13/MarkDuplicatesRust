/** Constructs a new FormatUtil and initializes various internal formatters.
 * This is necessary because SimpleDateFormat and other formatters are not threadsafe.
 */
pub(crate) struct FormatUtil {

}

impl FormatUtil {
    pub(crate) const DECIMAL_DIGITS_TO_PRINT:i32 = 6;
    
}
