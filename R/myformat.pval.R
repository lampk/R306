#' @export
myformat.pval <- function(p, cutoff = 0.0001){
  ## to format p value:
  ##- if <cutoff: report <cutoff; otherwise use only 2 decimal digits
  ## get number of digits: only work if the number of decimal places <15 (add 1 to avoid R mis-understanding)
  digits <- nchar(unlist(strsplit(as.character(cutoff + 1), split = "[.]"))[2])
  ## format
  out <- sapply(p, function(x) {
    if (is.na(x)) {
      NA
    } else {
      if (x >= cutoff) {
        format.pval(x, digits = 2, scientific = FALSE)
      } else {
        paste("<", formatC(cutoff, digits, format = "f"), sep = "")
      }
    }
  })
  return(out)
}