# Collects all functions from a directory
sourceDir <- function(path, trace = TRUE, ...) {
  op <- base::options(); on.exit(base::options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}