removeDepends <- function(pkg, recursive = FALSE){
  d <- package_dependencies(,installed.packages(), recursive = recursive)
  depends <- if(!is.null(d[[pkg]])) d[[pkg]] else character()
  needed <- unique(unlist(d[!names(d) %in% c(pkg,depends)]))
  toRemove <- depends[!depends %in% needed]
  if(length(toRemove)){
    toRemove <- select.list(c(pkg,sort(toRemove)), multiple = TRUE,
                            title = "Select packages to remove")
    remove.packages(toRemove)
    return(toRemove)
  } else {
    invisible(character())
  }
}
