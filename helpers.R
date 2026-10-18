## helpers.R -- two small functions the CAR model needs, copied from the
## wider project so this folder runs on its own.
##   season_of()        season label from a date ("2025-2026", July start)
##   write_inla_graph() writes the spatial adjacency for INLA's BYM2 term

## ---------------------------------------------------------------------------
## season helper (1 July boundary)
## ---------------------------------------------------------------------------
season_of <- function(week) {
  y <- lubridate::year(week); m <- lubridate::month(week)
  ifelse(m >= 7, paste0(y, "-", y + 1), paste0(y - 1, "-", y))
}

## ---------------------------------------------------------------------------
## adjacency -> INLA graph file
## ---------------------------------------------------------------------------
write_inla_graph <- function(W, file = tempfile(fileext = ".graph")) {
  n <- nrow(W)
  con <- file(file, "w"); on.exit(close(con))
  writeLines(as.character(n), con)
  for (i in seq_len(n)) {
    nb <- which(W[i, ] > 0)
    writeLines(paste(c(i, length(nb), nb), collapse = " "), con)
  }
  file
}
