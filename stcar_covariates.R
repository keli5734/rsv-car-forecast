
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(purrr); library(lubridate)
})

CACHE_DIR <- "data_cov"
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. POPULATION  (tidycensus -> state, population)
## ---------------------------------------------------------------------------
#' @param year  PEP vintage year
#' @param file  cache path
#' @param acs_fallback if the PEP product is unavailable for `year`, fall back to
#'   ACS 5-year total population (table B01003_001)
fetch_population <- function(year = 2023,
                             file = file.path(CACHE_DIR, "state_population.csv"),
                             acs_fallback = TRUE, overwrite = FALSE) {
  if (file.exists(file) && !overwrite) {
    message("population: using cached ", file)
    return(readr::read_csv(file, show_col_types = FALSE))
  }
  if (!requireNamespace("tidycensus", quietly = TRUE))
    stop("install.packages('tidycensus') and set a Census API key:\n",
         "  tidycensus::census_api_key('YOURKEY', install = TRUE)")

  pop <- try(
    tidycensus::get_estimates(geography = "state", product = "population",
                              vintage = year, year = year),
    silent = TRUE)

  if (inherits(pop, "try-error") || !nrow(as.data.frame(pop))) {
    if (!acs_fallback) stop("PEP request failed and acs_fallback = FALSE")
    message("population: PEP unavailable for ", year, " -> ACS 5-year B01003_001")
    pop <- tidycensus::get_acs(geography = "state", variables = "B01003_001",
                               year = year, survey = "acs5")
    out <- pop |> transmute(state = NAME, population = estimate)
  } else {
    pop <- as.data.frame(pop)
    vcol <- intersect(c("variable", "PRODUCT"), names(pop))
    if (length(vcol) && "POP" %in% pop[[vcol[1]]])
      pop <- pop[pop[[vcol[1]]] == "POP", ]
    out <- pop |> transmute(state = NAME, population = value)
  }

  out <- out |>
    filter(!state %in% c("Puerto Rico")) |>
    distinct(state, .keep_all = TRUE) |>
    arrange(state)
  readr::write_csv(out, file)
  message("population: wrote ", file, " (", nrow(out), " states)")
  out
}

## ---------------------------------------------------------------------------
## 2. CLIMATE  (Open-Meteo ERA5 archive -> state, week, temp, abs_humidity)
## ---------------------------------------------------------------------------
#' Absolute humidity (g/m^3) from temperature and dew point in degrees C.
#'   e  = 6.112 * exp(17.67 * Td / (Td + 243.5))      saturation vapour pressure
#'   AH = 2.1674 * e / (273.15 + T)
absolute_humidity <- function(temp_c, dew_c) {
  e <- 6.112 * exp(17.67 * dew_c / (dew_c + 243.5))
  2.1674 * e / (273.15 + temp_c)
}

#' Fetch daily weather at each state centroid and aggregate to ISO weeks.
#'
#' @param centroids data.frame(state, lat, lon); if NULL, taken from the bundled
#'   GeoJSON via state_centroids() in stcar_prep.R
#' @param start,end analysis window
#' @param pause seconds between API calls (be polite; the API is free)
fetch_climate <- function(centroids = NULL,
                          start = "2018-07-01", end = NULL,
                          file = file.path(CACHE_DIR, "state_climate.csv"),
                          pause = 1.5, max_retries = 5, resume = TRUE,
                          overwrite = FALSE) {

  ## Open-Meteo is free but rate-limited: 51 states x 8 years of daily data will
  ## trigger HTTP 429 if requested back-to-back. We therefore (a) pause between
  ## calls, (b) retry with exponential backoff on 429/5xx, and (c) SAVE
  ## INCREMENTALLY so an interrupted run can be resumed instead of restarted.

  if (file.exists(file) && !overwrite && !resume) {
    message("climate: using cached ", file)
    return(readr::read_csv(file, show_col_types = FALSE) |>
             mutate(week = as.Date(week)))
  }
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("install.packages('jsonlite')")

  if (is.null(centroids)) {
    if (exists("state_centroids")) {
      centroids <- try(state_centroids(), silent = TRUE)
    } else centroids <- NULL
    if (is.null(centroids) || inherits(centroids, "try-error")) {
      message("climate: using built-in state centres (no geojson available)")
      centroids <- data.frame(state = datasets::state.name,
                              lon = datasets::state.center$x,
                              lat = datasets::state.center$y,
                              stringsAsFactors = FALSE)
      centroids$lon[centroids$state == "Alaska"] <- -152.4
      centroids$lat[centroids$state == "Alaska"] <- 64.0
      centroids$lon[centroids$state == "Hawaii"] <- -155.5
      centroids$lat[centroids$state == "Hawaii"] <- 19.6
      centroids <- rbind(centroids,
        data.frame(state = "District of Columbia", lon = -77.03, lat = 38.90))
    }
  }
  stopifnot(all(c("state", "lat", "lon") %in% names(centroids)))

  ## the ERA5 archive lags real time; never ask beyond a few days back
  max_end <- as.character(Sys.Date() - 7)
  if (is.null(end) || as.Date(end) > as.Date(max_end)) end <- max_end

  ## ---- resume: keep whatever already succeeded -----------------------------
  have <- NULL
  if (resume && file.exists(file) && !overwrite) {
    have <- readr::read_csv(file, show_col_types = FALSE) |> mutate(week = as.Date(week))
    done <- unique(have$state)
    todo <- setdiff(centroids$state, done)
    message("climate: resuming - ", length(done), " states cached, ",
            length(todo), " to fetch")
  } else {
    todo <- centroids$state
  }
  if (!length(todo)) {
    message("climate: all states already cached")
    return(have)
  }

  base <- "https://archive-api.open-meteo.com/v1/archive"
  daily_vars <- "temperature_2m_mean,dew_point_2m_mean,precipitation_sum"

  get_one <- function(st) {
    i <- match(st, centroids$state)
    url <- sprintf(
      "%s?latitude=%.4f&longitude=%.4f&start_date=%s&end_date=%s&daily=%s&timezone=UTC",
      base, centroids$lat[i], centroids$lon[i], start, end, daily_vars)
    for (attempt in seq_len(max_retries)) {
      js <- try(suppressWarnings(jsonlite::fromJSON(url)), silent = TRUE)
      if (!inherits(js, "try-error") && !is.null(js$daily)) {
        return(tibble::tibble(state = st,
                              date   = as.Date(js$daily$time),
                              temp   = js$daily$temperature_2m_mean,
                              dew    = js$daily$dew_point_2m_mean,
                              precip = js$daily$precipitation_sum))
      }
      wait <- pause * 2^(attempt - 1) + runif(1, 0, 1)   # exponential backoff
      message("  ", st, ": attempt ", attempt, " failed (likely HTTP 429) - ",
              "waiting ", round(wait, 1), "s")
      Sys.sleep(wait)
    }
    warning("climate: giving up on ", st,
            " after ", max_retries, " attempts. Re-run fetch_climate() later; ",
            "progress is cached and the run will resume.")
    NULL
  }

  message("climate: querying Open-Meteo for ", length(todo), " states (",
          start, " .. ", end, "); pause = ", pause, "s, retries = ", max_retries)

  acc <- have
  for (k in seq_along(todo)) {
    st <- todo[k]
    raw1 <- get_one(st)
    Sys.sleep(pause)
    if (is.null(raw1)) next
    wk1 <- raw1 |>
      mutate(abs_humidity = absolute_humidity(temp, dew),
             week = as.Date(cut(date, "week", start.on.monday = TRUE)) + 5) |>
      group_by(state, week) |>
      summarise(temp         = mean(temp, na.rm = TRUE),
                dew_point    = mean(dew,  na.rm = TRUE),
                abs_humidity = mean(abs_humidity, na.rm = TRUE),
                precip       = sum(precip, na.rm = TRUE),
                .groups = "drop")
    acc <- dplyr::bind_rows(acc, wk1)
    readr::write_csv(acc, file)          # save after EVERY state
    if (k %% 10 == 0 || k == length(todo))
      message("  ", k, "/", length(todo), " states done (cached)")
  }

  if (is.null(acc) || !nrow(acc)) stop("climate: no data returned")
  missing <- setdiff(centroids$state, unique(acc$state))
  if (length(missing))
    message("climate: STILL MISSING ", length(missing), " states (",
            paste(head(missing, 5), collapse = ", "),
            if (length(missing) > 5) ", ..." else "",
            ") - re-run fetch_climate() to resume")
  message("climate: wrote ", file, " (", nrow(acc), " state-weeks, ",
          dplyr::n_distinct(acc$state), " states)")
  acc
}

## ---------------------------------------------------------------------------
## 3. MOBILITY  (ACS state-to-state migration flows -> origin, destination, flow)
## ---------------------------------------------------------------------------
#' @param year ACS 5-year endyear; the Migration Flows API covers 2010-2018
fetch_mobility <- function(year = 2018,
                           file = file.path(CACHE_DIR, "state_flows.csv"),
                           states = NULL, overwrite = FALSE, pause = 0.2) {
  if (file.exists(file) && !overwrite) {
    message("mobility: using cached ", file)
    return(readr::read_csv(file, show_col_types = FALSE))
  }
  if (!requireNamespace("tidycensus", quietly = TRUE))
    stop("install.packages('tidycensus') and set a Census API key")

  ## NOTE: the ACS Migration Flows API does NOT support geography = "state".
  ## It serves "county", "county subdivision" and "metropolitan statistical
  ## area" only. We therefore request COUNTY flows one state at a time and
  ## aggregate them up to a state x state matrix. That is ~51 API calls, so the
  ## result is cached.
  if (is.null(states)) states <- c(datasets::state.name, "District of Columbia")

  message("mobility: ACS county flows for ", length(states),
          " states (", year, "), aggregating to state x state ...")

  pull_one <- function(st) {
    fl <- try(tidycensus::get_flows(geography = "county", state = st,
                                    year = year, output = "wide"),
              silent = TRUE)
    if (inherits(fl, "try-error")) {
      warning("mobility: failed for ", st, " - ", conditionMessage(attr(fl, "condition")))
      return(NULL)
    }
    Sys.sleep(pause)
    nm <- names(fl)
    ## FULL2_NAME is the OTHER end of the flow, e.g. "Cook County, Illinois"
    d_col <- intersect(c("FULL2_NAME", "NAME2"), nm)[1]
    v_col <- intersect(c("MOVEDIN", "MOVEDNET", "MOVEDOUT"), nm)[1]
    if (is.na(d_col) || is.na(v_col)) return(NULL)
    tibble::tibble(destination = st,
                   other = as.character(fl[[d_col]]),
                   flow  = suppressWarnings(as.numeric(fl[[v_col]])))
  }

  raw <- purrr::map_dfr(states, pull_one)
  if (!nrow(raw)) stop("mobility: no flows returned (check your Census API key)")

  ## "County, State" -> "State"; drop foreign origins (no comma / not a state)
  valid <- c(datasets::state.name, "District of Columbia")
  out <- raw |>
    mutate(origin = trimws(sub(".*,\\s*", "", other))) |>
    filter(!is.na(flow), flow > 0,
           origin %in% valid, destination %in% valid,
           origin != destination) |>
    group_by(origin, destination) |>
    summarise(flow = sum(flow, na.rm = TRUE), .groups = "drop")

  ## symmetrise: connectivity should not depend on migration direction
  out <- out |>
    bind_rows(out |> rename(origin = destination, destination = origin)) |>
    group_by(origin, destination) |>
    summarise(flow = sum(flow), .groups = "drop")

  readr::write_csv(out, file)
  message("mobility: wrote ", file, " (", nrow(out), " state pairs)")
  out
}

## ---------------------------------------------------------------------------
## 3b. COMMUTING FLOWS  (LODES -> origin, destination, flow)
## ---------------------------------------------------------------------------
#' State-to-state COMMUTING flows from LEHD LODES.
#'
#' WHY THIS AND NOT MIGRATION. `fetch_mobility()` returns ACS Migration Flows,
#' which count people who CHANGED RESIDENCE -- permanent moves, averaged over
#' 2014-2018. Weekly respiratory transmission is driven by people who travel
#' repeatedly between two places, which is journey-to-work, not migration.
#'
#' LODES records, for every job, the census block of the WORKPLACE and the
#' block of the worker's RESIDENCE. The "aux" part of each state's origin-
#' destination file holds exactly the jobs whose workplace is in that state and
#' whose worker lives OUTSIDE it -- cross-state commuting, which is what we
#' want. Aggregating the home block's state FIPS gives an origin-by-destination
#' commuter matrix.
#'
#' Interpretation: flow[i, j] is the number of people who LIVE in state i and
#' WORK in state j. It is directional, and unlike migration it is a recurring
#' movement, so it is a far better proxy for regular contact between states.
#'
#' LIMITATIONS to state in the paper:
#'   * jobs only -- it misses non-work travel (family visits, holidays), which
#'     matters for a paediatric infection;
#'   * one vintage, so the matrix is time-invariant, as migration was;
#'   * LODES excludes most federal employment and some self-employment.
#'
#' @param year   LODES year; 2021 is widely available in LODES8
#' @param version LODES version passed to lehdr
fetch_commuting <- function(year = 2021,
                            file = file.path(CACHE_DIR, "state_commute.csv"),
                            states = NULL, overwrite = FALSE,
                            lodes_version = "LODES8") {

  if (file.exists(file) && !overwrite) {
    message("commuting: using cached ", file)
    return(readr::read_csv(file, show_col_types = FALSE))
  }
  if (!requireNamespace("lehdr", quietly = TRUE))
    stop("install.packages('lehdr')  -- it downloads LODES from the Census",
         call. = FALSE)

  if (is.null(states)) states <- c(datasets::state.name, "District of Columbia")
  abb <- c(setNames(datasets::state.abb, datasets::state.name),
           "District of Columbia" = "DC")
  st_abb <- tolower(unname(abb[states]))
  st_abb <- st_abb[!is.na(st_abb)]

  message("commuting: downloading LODES ", year, " cross-state OD for ",
          length(st_abb), " states (this takes a while on the first run)")

  ## FIPS -> state name, for mapping the home block back to a state
  fips2name <- setNames(names(FIPS_STATE), FIPS_STATE)      # name -> fips chr
  name_by_fips <- setNames(FIPS_STATE, sprintf("%02d", as.integer(names(FIPS_STATE))))

  pull_one <- function(ab) {
    ## state_part = "aux": workplace in this state, residence OUTSIDE it
    od <- try(lehdr::grab_lodes(state = ab, year = year,
                                lodes_type = "od", job_type = "JT00",
                                segment = "S000", state_part = "aux",
                                agg_geo = "block", version = lodes_version),
              silent = TRUE)
    if (inherits(od, "try-error")) {
      warning("commuting: failed for ", toupper(ab), " - ",
              conditionMessage(attr(od, "condition")), call. = FALSE)
      return(NULL)
    }
    hc <- grep("^h_geocode", names(od), value = TRUE)[1]
    jc <- if ("S000" %in% names(od)) "S000" else
            grep("^S000", names(od), value = TRUE)[1]
    if (is.na(hc) || is.na(jc)) return(NULL)
    tibble::tibble(
      dest_abb  = toupper(ab),
      home_fips = substr(as.character(od[[hc]]), 1, 2),
      jobs      = as.numeric(od[[jc]])) |>
      dplyr::group_by(dest_abb, home_fips) |>
      dplyr::summarise(jobs = sum(jobs, na.rm = TRUE), .groups = "drop")
  }

  raw <- purrr::map_dfr(st_abb, pull_one)
  if (!nrow(raw)) stop("commuting: no LODES data returned", call. = FALSE)

  name_by_abb <- setNames(names(abb), unname(abb))
  out <- raw |>
    mutate(origin      = unname(name_by_fips[home_fips]),
           destination = unname(name_by_abb[dest_abb])) |>
    filter(!is.na(origin), !is.na(destination), origin != destination) |>
    group_by(origin, destination) |>
    summarise(flow = sum(jobs, na.rm = TRUE), .groups = "drop")

  ## Commuting is DIRECTIONAL: living in NJ and working in NY is not the same
  ## trip as the reverse. For a symmetric connectivity matrix we add the two
  ## directions, which is the total number of people crossing that border for
  ## work in either direction.
  out <- out |>
    bind_rows(out |> rename(origin = destination, destination = origin)) |>
    group_by(origin, destination) |>
    summarise(flow = sum(flow), .groups = "drop")

  readr::write_csv(out, file)
  message("commuting: wrote ", file, " (", nrow(out), " state pairs, ",
          format(sum(out$flow), big.mark = ","), " commuters)")
  out
}


## ---------------------------------------------------------------------------
## convenience: fetch all three
## ---------------------------------------------------------------------------
fetch_all_covariates <- function(start = "2018-07-01", end = NULL,
                                 pop_year = 2023, flow_year = 2018,
                                 connectivity = c("commute", "migration"),
                                 lodes_year = 2021,
                                 centroids = NULL, overwrite = FALSE) {
  connectivity <- match.arg(connectivity)
  pop <- try(fetch_population(year = pop_year, overwrite = overwrite), silent = TRUE)
  if (inherits(pop, "try-error"))
    message("!! population unavailable: ", conditionMessage(attr(pop, "condition")))
  ## "commute" (LODES journey-to-work) is the right measure for recurring
  ## contact between states; "migration" (ACS) counts permanent moves and is
  ## kept only for comparison.
  mob <- if (connectivity == "commute")
    try(fetch_commuting(year = lodes_year, overwrite = overwrite), silent = TRUE)
  else
    try(fetch_mobility(year = flow_year, overwrite = overwrite), silent = TRUE)
  if (inherits(mob, "try-error"))
    message("!! mobility unavailable: ", conditionMessage(attr(mob, "condition")))
  clim <- try(fetch_climate(centroids = centroids, start = start, end = end,
                            overwrite = overwrite), silent = TRUE)
  if (inherits(clim, "try-error"))
    message("!! climate unavailable: ", conditionMessage(attr(clim, "condition")))

  list(
    pop_file      = if (!inherits(pop,  "try-error")) file.path(CACHE_DIR, "state_population.csv") else NULL,
    mobility_file = if (!inherits(mob, "try-error"))
      file.path(CACHE_DIR, if (connectivity == "commute") "state_commute.csv"
                           else "state_flows.csv") else NULL,
    climate_file  = if (!inherits(clim, "try-error")) file.path(CACHE_DIR, "state_climate.csv")    else NULL
  )
}
