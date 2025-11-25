### Create a table for the package manual

## https://roxygen2.r-lib.org/articles/formatting.html
tabular <- function(df, ...) {
  stopifnot(is.data.frame(df))

  align <- function(x) if (is.numeric(x)) "r" else "l"
  col_align <- vapply(df, align, character(1))

  cols <- lapply(df, format, ...)
  contents <- do.call(
    "paste",
    c(cols, list(sep = " \\tab ", collapse = "\\cr\n   "))
  )

  paste(
    "\\tabular{",
    paste(col_align, collapse = ""),
    "}{\n   ",
    paste0("\\strong{", names(df), "}", sep = "", collapse = " \\tab "),
    " \\cr\n   ",
    contents,
    "\n }\n",
    sep = ""
  )
}


tab_dists <- function() {
  ns <- asNamespace("energyGOF")
  dists <- lsf.str(envir = ns)
  dists <- dists[grepl("_dist$", dists)]
  dists <- dists[!(dists %in% c("xform_dist", "char_to_dist"))]

  d <- do.call(dists[1], args = as.list(formals(dists[1])))
  names(d$par)
  d$composite_p

  dat <- data.frame(
    Distribution = character(),
    Function = character(),
    Parameters = character(),
    Composite_Test = logical()
  )

  for (i in seq_along(dists)) {
    d <- do.call(dists[i], args = as.list(formals(dists[i])))
    Distribution = d$name
    Function = dists[i]
    Parameters = paste(names(d$par), collapse = ", ")
    Composite_Test = d$composite_p
    row <- data.frame(
      Distribution = Distribution,
      Function = Function,
      Parameters = Parameters,
      Composite_Test = Composite_Test
    )
    dat <- rbind(dat, row)
  }
  tabular(dat)
}
