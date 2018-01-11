#' @export
myformat.pval <- function(p, cutoff = 0.001, type = 1){
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
        if (type == 1) {
          formatC(x, digits, format = "f")
        } else {
          if (type == 2) {
            format.pval(x, digits = 2, scientific = FALSE)
          }
        }
        
      } else {
        paste("<", formatC(cutoff, digits, format = "f"), sep = "")
      }
    }
  })
  return(out)
}