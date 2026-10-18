## =============================================================================
## hub_data.R — the CDC RSV Forecast Hub target data, in this project's shape
##
##   source("R/hub_data.R")
##   hub <- hub_target_data()        # downloads and reshapes
##
## ---------------------------------------------------------------------------
## THE TARGET IS NOT EPIC COSMOS
## ---------------------------------------------------------------------------
## This matters more than any formatting detail. The hub collects forecasts of
##
##   wk inc rsv hosp             NHSN laboratory-confirmed RSV hospital
##                               ADMISSIONS (mandatory)
##   wk inc rsv prop ed visits   NSSP proportion of ED visits due to RSV,
##                               as a decimal (optional)
##
## Epic Cosmos is neither. Two consequences:
##
##   1. The models must be REFITTED on the hub's series. Nothing estimated from
##      Epic transfers -- not the state effects, not the seasonal amplitudes,
##      not rho.
##
##   2. NEITHER TARGET IS CENSORED. The interval-censored likelihood that makes
##      the Bayesian model distinctive here has nothing to do in a hub
##      submission. The CAR structure, the harmonics and the AR(1) all still
##      apply; the censoring simply switches off. That is worth saying plainly
##      rather than discovering at the point of writing the methods.
##
## What DOES transfer is the comparison itself: three structurally different
## models and a combination rule, tested against each other on identical
## targets. That is the contribution, and it is portable.
##
## ---------------------------------------------------------------------------
## HORIZONS RUN FROM -1
## ---------------------------------------------------------------------------
##   -1  the epiweek BEFORE the reference date -- already past, but not yet
##       fully reported
##    0  the current epiweek
##    1, 2, 3  the next three
##
## So two of the five are nowcasts of weeks that have happened. That is the
## reporting-delay problem, and it is where this project's work is strongest.
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

HUB_RAW <- "https://raw.githubusercontent.com/CDCgov/rsv-forecast-hub/main"

#' Download the hub's target data and location table.
#'
#' @param target  "hosp" for NHSN admissions (the mandatory target) or "ed" for
#'   the NSSP ED proportion
#' @param cache   where to keep the downloads, so a rerun does not refetch
hub_target_data <- function(target = c("hosp", "ed"), cache = "data_hub",
                            refresh = FALSE, verbose = TRUE) {
  target <- match.arg(target)
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)

  ## ---- locations ---------------------------------------------------------
  lf <- file.path(cache, "locations.csv")
  if (refresh || !file.exists(lf))
    download.file(file.path(HUB_RAW, "auxiliary-data/locations.csv"), lf,
                  quiet = !verbose)
  loc <- read_csv(lf, show_col_types = FALSE)
  if (verbose) {
    cat("\n=========== LOCATIONS ===========\n")
    cat(nrow(loc), " rows;  columns: ", paste(names(loc), collapse = ", "),
        "\n", sep = "")
    print(head(as.data.frame(loc), 4), row.names = FALSE)
  }

  ## ---- target data -------------------------------------------------------
  ## PARQUET, not CSV: the hub publishes target-data/time-series.parquet, which
  ## holds BOTH targets in long form. `arrow` is therefore required.
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("install.packages('arrow') -- the hub's target data is parquet",
         call. = FALSE)

  tf <- file.path(cache, "time-series.parquet")
  if (refresh || !file.exists(tf))
    download.file(file.path(HUB_RAW, "target-data/time-series.parquet"), tf,
                  mode = "wb", quiet = !verbose)

  td <- arrow::read_parquet(tf)
  if (verbose) {
    cat("\n=========== TARGET DATA ===========\n")
    cat("columns: ", paste(names(td), collapse = ", "), "\n", sep = "")
    print(head(as.data.frame(td), 4), row.names = FALSE)
  }

  ## ---- reshape to the project's panel ------------------------------------
  ## columns: target_end_date, as_of, location, target, observation
  tname <- if (target == "hosp") "wk inc rsv hosp" else
    "wk inc rsv prop ed visits"

  ## AS_OF IS A VINTAGE. The hub ships the full reporting triangle: each
  ## (week, location) appears once per data release, so the file records what
  ## was known at every past moment. Two consequences.
  ##
  ## For FITTING, only the latest vintage should be used, or the same week
  ## enters many times over.
  ##
  ## For BACKTESTING, this is what makes a real-time replay possible: a
  ## forecast made at a past reference date can be given exactly the data that
  ## existed then, rather than today's revised figures. Scoring against
  ## revised data while fitting on revised data flatters any model, and the
  ## hub's own evaluation does not work that way.
  ## `hub_vintage()` below extracts a single past snapshot.
  d <- td |>
    filter(target == tname) |>
    mutate(week = as.Date(target_end_date),
           as_of = as.Date(as_of),
           location = as.character(location),
           y = suppressWarnings(as.numeric(observation)))

  if (verbose)
    cat("\nvintages: ", n_distinct(d$as_of), " releases, ",
        format(min(d$as_of)), " to ", format(max(d$as_of)),
        "\n  (keeping the latest for fitting; see hub_vintage() to replay a ",
        "past one)\n", sep = "")
  d <- d |> filter(as_of == max(as_of))

  ## state names, so the existing code can be reused unchanged
  lnm <- intersect(c("location_name", "abbreviation", "state_name"), names(loc))[1]
  d <- d |>
    left_join(loc |> transmute(location = as.character(location),
                               state = .data[[lnm]]),
              by = "location") |>
    filter(!is.na(week), is.finite(y)) |>
    select(state, location, week, y) |>
    arrange(state, week)

  if (verbose) {
    cat("\npanel: ", nrow(d), " rows;  ", n_distinct(d$state), " locations;  ",
        format(min(d$week)), " to ", format(max(d$week)),
        "  (", n_distinct(d$week), " weeks, about ",
        round(n_distinct(d$week) / 52, 1), " seasons)\n", sep = "")
    if (target == "hosp")
      cat("  NHSN began collecting RSV admissions in late 2023, so this is ",
          "the whole\n  series. The ED-proportion target reaches back to ",
          "2022-10 -- a year more --\n  and is the better one to develop on.",
          "\n", sep = "")
    cat("value: median ", sprintf("%.2f", median(d$y)),
        ", max ", sprintf("%.0f", max(d$y)), "\n", sep = "")
    cat("\nNOTE: this series is NOT censored. The models are refitted here ",
        "with the\ncensored likelihood switched off.\n", sep = "")
  }
  list(data = d, locations = loc, target = target, file = tf)
}


#' The data as it stood at a past release.
#'
#' A backtest that fits on today's revised figures and scores against today's
#' revised figures gives every model information it could not have had. The
#' hub's own evaluation compares a forecast made on a Wednesday against what
#' was eventually reported, having been made from what was known that day.
#' This reproduces that: pass a past `as_of` and get the file exactly as it
#' was then.
#'
#' NHSN admissions revise little -- the first report is about 93% of the final
#' value one week out and complete by two -- so the effect is small for that
#' target. It is worth checking rather than assuming, and it matters more for
#' the ED proportion.
hub_vintage <- function(cache = "data_hub", as_of,
                        target = c("hosp", "ed"), locations = NULL) {
  target <- match.arg(target)
  tf <- file.path(cache, "time-series.parquet")
  if (!file.exists(tf)) stop("run hub_target_data() first", call. = FALSE)
  tname <- if (target == "hosp") "wk inc rsv hosp" else
    "wk inc rsv prop ed visits"
  d <- arrow::read_parquet(tf) |>
    filter(target == tname, as.Date(as_of) <= as.Date(!!as_of)) |>
    mutate(week = as.Date(target_end_date), location = as.character(location),
           y = as.numeric(observation))
  if (!is.null(locations)) d <- filter(d, location %in% locations)
  ## the most recent figure for each week AS OF that date
  d |> group_by(location, week) |>
    slice_max(as.Date(as_of), n = 1, with_ties = FALSE) |>
    ungroup() |> select(location, week, y) |> arrange(location, week)
}


#' How much does the target revise after first publication?
#'
#' Worth running once per target. If the first report is already close to the
#' final value there is nothing for a delay model to correct, and a nowcast at
#' horizon -1 reduces to reading the data. If it is not, horizons -1 and 0 are
#' where a model can add the most.
hub_revision <- function(cache = "data_hub", target = c("hosp", "ed"),
                         max_delay = 8, min_final = 5) {
  target <- match.arg(target)
  tf <- file.path(cache, "time-series.parquet")
  tname <- if (target == "hosp") "wk inc rsv hosp" else
    "wk inc rsv prop ed visits"
  d <- arrow::read_parquet(tf) |>
    filter(target == tname) |>
    mutate(week = as.Date(target_end_date), as_of = as.Date(as_of),
           location = as.character(location), y = as.numeric(observation)) |>
    mutate(delay = as.integer(round(as.numeric(as_of - week) / 7))) |>
    filter(delay >= 0, delay <= max_delay)

  fin <- d |> filter(as_of == max(as_of), y > min_final) |>
    select(location, week, final = y)
  cur <- d |> inner_join(fin, by = c("location", "week")) |>
    mutate(p = y / final)

  out <- cur |> group_by(delay) |>
    summarise(n = n(), completeness = median(p), .groups = "drop")
  cat("\n=========== REVISION IN THE HUB TARGET (", target,
      ") ===========\n", sep = "")
  print(as.data.frame(out), row.names = FALSE, digits = 3)
  cat("\nIf completeness reaches 1 within a week or two, horizon -1 is ",
      "nearly an\nobservation and the nowcast adds little there. The longer ",
      "the tail, the\nmore a delay-aware model can contribute.\n", sep = "")
  invisible(out)
}


#' Put the hub panel into the shape `stcar_fit_inla3()` expects.
#'
#' THE OFFSET IS POPULATION. The Epic panel carried ED encounters as the
#' denominator; NHSN reports admissions with no denominator at all, so the
#' natural exposure is the state's population, which the hub ships in
#' `locations.csv`. Without it the state effect b_i would have to absorb the
#' whole size difference between California and Wyoming -- it can, but then
#' b_i measures population rather than anything about RSV, and the spatial
#' prior is smoothing the wrong quantity between neighbours.
#'
#' With population as the offset the model is on admissions per person, which
#' is the scale on which "neighbouring states are similar" is a reasonable
#' thing to assume.
#'
#' For the ED-proportion target the observation is already a proportion. It is
#' scaled to a pseudo-count against a fixed denominator so the count models
#' apply; the Poisson variance is then an assumption rather than a property of
#' the data, which is worth stating in the methods.
hub_to_panel <- function(hub, drop_locations = c("US", "72"),
                         ed_scale = 1e5, per = 1e5, verbose = TRUE) {
  d <- hub$data |> filter(!location %in% drop_locations)
  states <- sort(unique(d$state))

  pop <- hub$locations |>
    transmute(location = as.character(location),
              pop = as.numeric(population))
  d <- d |> left_join(pop, by = "location")
  if (anyNA(d$pop))
    stop("no population for: ",
         paste(unique(d$location[is.na(d$pop)]), collapse = ", "),
         call. = FALSE)

  if (hub$target == "ed") {
    d <- d |> mutate(encounters = ed_scale, y_raw = round(y * ed_scale))
    if (verbose)
      message("ED proportions scaled by ", format(ed_scale, big.mark = ","),
              " to pseudo-counts; the Poisson variance is then an assumption, ",
              "not a property of the data")
  } else {
    ## population in units of `per`, so the intercept is log admissions per
    ## 100,000 and is of a readable magnitude
    d <- d |> mutate(encounters = pop / per * 1e5, y_raw = round(y))
    if (verbose)
      cat("offset: state population (per ", format(per, big.mark = ","),
          ");  range ", format(round(min(d$pop) / 1e6, 2)), "M to ",
          format(round(max(d$pop) / 1e6, 1)), "M\n", sep = "")
  }

  out <- d |>
    transmute(state, week, y_raw, encounters,
              censored = FALSE,          # the hub targets are not suppressed
              log_offset = log(encounters))
  if (verbose)
    cat("\npanel for fitting: ", nrow(out), " rows, ", length(states),
        " locations, censored = FALSE throughout\n", sep = "")
  list(data = out, states = states, locations = unique(d$location))
}
