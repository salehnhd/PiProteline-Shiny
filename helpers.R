# helpers.R -------------------------------------------------------------------
# Small utilities shared by every tab. Sourced once from app.R BEFORE the tabs,
# so these names are available to all server functions at run time.

# Null/empty coalescing operator. The old code used `%||%` without defining or
# importing it, which can fail (shiny does not export it). Define it here once.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Safe nested extraction: get_path(list, c("a","b","c")) -> list$a$b$c or NULL.
get_path <- function(x, path) {
  out <- x
  for (nm in path) {
    if (is.null(out) || is.null(out[[nm]])) return(NULL)
    out <- out[[nm]]
  }
  out
}

# Build human-readable labels for a list of pairwise results (volcano / MDS).
# Uses the element names if the backend provided them; otherwise falls back to
# the same group order combn() would produce, so labels match any group count.
contrast_labels <- function(lst, groups) {
  if (is.null(lst) || length(lst) == 0) return(character(0))
  nm <- names(lst)
  if (!is.null(nm) && all(nzchar(nm))) return(nm)
  if (length(groups) >= 2) {
    combos <- utils::combn(groups, 2, function(x) paste(x[1], "vs", x[2]))
    if (length(combos) == length(lst)) return(combos)
  }
  paste0("contrast_", seq_along(lst))
}
