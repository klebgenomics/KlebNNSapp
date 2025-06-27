# This helper function filters the dots and returns only those arguments that
# belong to the internal function.

match.args <- function(dots, fn) {
  fn_args <- names(formals(fn))
  matched_args <- dots[names(dots) %in% fn_args]
  return(matched_args)
}

exec2 <- name <- function(fn, ..., args=list()) {
  list(...) |>
    match.args(fn) |>
    utils::modifyList(args) |>
    do.call(what=fn)
}
