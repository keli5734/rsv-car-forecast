
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(stringr); library(lubridate)
})

if (file.exists("R/DataExtract.R")) source("R/DataExtract.R") else
if (file.exists("DataExtract.R"))   source("DataExtract.R")   else
if (file.exists("R/data.R"))        source("R/data.R")        else
if (file.exists("data.R"))          source("data.R")
## The Epic Cosmos data layer is not needed for the NHSN CAR model, which
## builds its panel with stcar_panel(); season_of() comes from helpers.R.

## real covariate fetchers (population / climate / mobility)
if (file.exists("R/stcar_covariates.R")) source("R/stcar_covariates.R") else
if (file.exists("stcar_covariates.R"))   source("stcar_covariates.R")

CENSOR_LO <- 1L    # a suppressed week means the true count is in [CENSOR_LO,
CENSOR_HI <- 10L   # CENSOR_HI]; genuine counts begin at 11

## ---------------------------------------------------------------------------
## state adjacency (shared land border) -> W for the CAR prior
## ---------------------------------------------------------------------------
ADJ_USPS <- c(
  AL="FL GA MS TN", AZ="CA CO NM NV UT", AR="LA MO MS OK TN TX", CA="AZ NV OR",
  CO="AZ KS NE NM OK UT WY", CT="MA NY RI", DE="MD NJ PA", DC="MD VA",
  FL="AL GA", GA="AL FL NC SC TN", ID="MT NV OR UT WA WY", IL="IN IA KY MO WI",
  IN="IL KY MI OH", IA="IL MN MO NE SD WI", KS="CO MO NE OK",
  KY="IL IN MO OH TN VA WV", LA="AR MS TX", ME="NH", MD="DC DE PA VA WV",
  MA="CT NH NY RI VT", MI="IN OH WI", MN="IA ND SD WI", MS="AL AR LA TN",
  MO="AR IA IL KS KY NE OK TN", MT="ID ND SD WY", NE="CO IA KS MO SD WY",
  NV="AZ CA ID OR UT", NH="MA ME VT", NJ="DE NY PA", NM="AZ CO OK TX UT",
  NY="CT MA NJ PA VT", NC="GA SC TN VA", ND="MN MT SD", OH="IN KY MI PA WV",
  OK="AR CO KS MO NM TX", OR="CA ID NV WA", PA="DE MD NJ NY OH WV",
  RI="CT MA", SC="GA NC", SD="IA MN MT ND NE WY",
  TN="AL AR GA KY MO MS NC VA", TX="AR LA NM OK", UT="AZ CO ID NM NV WY",
  VT="MA NH NY", VA="DC KY MD NC TN WV", WA="ID OR", WV="KY MD OH PA VA",
  WI="IA IL MI MN", WY="CO ID MT NE SD UT")

USPS_NAME <- c(
  AL="Alabama", AK="Alaska", AZ="Arizona", AR="Arkansas", CA="California",
  CO="Colorado", CT="Connecticut", DE="Delaware", DC="District of Columbia",
  FL="Florida", GA="Georgia", HI="Hawaii", ID="Idaho", IL="Illinois",
  IN="Indiana", IA="Iowa", KS="Kansas", KY="Kentucky", LA="Louisiana",
  ME="Maine", MD="Maryland", MA="Massachusetts", MI="Michigan", MN="Minnesota",
  MS="Mississippi", MO="Missouri", MT="Montana", NE="Nebraska", NV="Nevada",
  NH="New Hampshire", NJ="New Jersey", NM="New Mexico", NY="New York",
  NC="North Carolina", ND="North Dakota", OH="Ohio", OK="Oklahoma", OR="Oregon",
  PA="Pennsylvania", RI="Rhode Island", SC="South Carolina", SD="South Dakota",
  TN="Tennessee", TX="Texas", UT="Utah", VT="Vermont", VA="Virginia",
  WA="Washington", WV="West Virginia", WI="Wisconsin", WY="Wyoming")

#' Binary adjacency matrix over a given ordered vector of state names.
#' Alaska and Hawaii have no land neighbours; they are connected to their
#' nearest mainland state so the CAR graph stays connected (islands otherwise
#' break the intrinsic CAR prior).
build_adjacency <- function(states, connect_islands = TRUE) {
  n <- length(states); W <- matrix(0L, n, n, dimnames = list(states, states))
  name2usps <- setNames(names(USPS_NAME), USPS_NAME)
  for (s in states) {
    u <- name2usps[[s]]
    if (is.null(u) || is.na(u) || !u %in% names(ADJ_USPS)) next
    nb <- unlist(strsplit(ADJ_USPS[[u]], " "))
    nb <- unname(USPS_NAME[nb]); nb <- nb[nb %in% states]
    if (length(nb)) W[s, nb] <- 1L
  }
  W <- pmax(W, t(W))                                   # force symmetry
  if (connect_islands) {
    if (all(c("Alaska", "Washington") %in% states))
      W["Alaska", "Washington"] <- W["Washington", "Alaska"] <- 1L
    if (all(c("Hawaii", "California") %in% states))
      W["Hawaii", "California"] <- W["California", "Hawaii"] <- 1L
  }
  stopifnot(all(rowSums(W) > 0))                       # no isolated nodes
  W
}

## ---------------------------------------------------------------------------
## DISTANCE-BASED spatial weights (alternative to binary adjacency)
## ---------------------------------------------------------------------------
## Binary adjacency has three weaknesses for states:
##   * it treats New Hampshire -> Maine identically to Texas -> Oklahoma;
##   * state size is arbitrary relative to epidemiology, so a large state has
##     the same single link to a neighbour as a small one;
##   * Alaska and Hawaii have no land border at all and must be attached by hand.
##
## A distance-decay graph on state centroids avoids all three. The ICAR prior in
## NIMBLE accepts general non-negative symmetric weights, so
##     u_i | u_-i ~ N( sum_j w_ij u_j / sum_j w_ij ,  1 / sum_j w_ij )
## replaces the unweighted mean over land neighbours.
##
## The graph must stay SPARSE. A fully dense distance matrix connects every
## state to every other and the ICAR degenerates into a global mean, removing
## the spatial structure it is meant to encode. We therefore keep either the k
## nearest neighbours or a distance cut-off, then apply the decay.
##
## Decay options (d in km):
##   "exp"      w = exp(-d / rho)              rho = e-folding distance
##   "inverse"  w = 1 / (d / rho)              heavier tail than exponential
##   "gaussian" w = exp(-(d / rho)^2)          sharper local falloff
##
build_distance_W <- function(states, centroids,
                             type = c("knn", "cutoff"),
                             k = 6, cutoff_km = 1200,
                             decay = c("exp", "inverse", "gaussian"),
                             rho_km = 600, standardise = TRUE) {
  type  <- match.arg(type)
  decay <- match.arg(decay)
  cc <- centroids[match(states, centroids$state), ]
  if (anyNA(cc$lat))
    stop("missing centroid for: ", paste(states[is.na(cc$lat)], collapse = ", "),
         "\n  state_centroids() fills District of Columbia, Alaska and Hawaii",
         "\n  automatically; any other gap means the GeoJSON uses different",
         "\n  state names. Pass geojson = NULL to use the built-in centres.",
         call. = FALSE)
  n <- length(states)

  ## great-circle distance (km)
  rad <- function(x) x * pi / 180
  D <- matrix(0, n, n, dimnames = list(states, states))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j) next
    dlat <- rad(cc$lat[j] - cc$lat[i]); dlon <- rad(cc$lon[j] - cc$lon[i])
    a <- sin(dlat/2)^2 + cos(rad(cc$lat[i])) * cos(rad(cc$lat[j])) * sin(dlon/2)^2
    D[i, j] <- 2 * 6371 * asin(pmin(1, sqrt(a)))
  }

  ## sparsify
  M <- matrix(0, n, n, dimnames = list(states, states))
  if (type == "knn") {
    for (i in seq_len(n)) {
      o <- order(D[i, ]); o <- o[o != i][seq_len(min(k, n - 1))]
      M[i, o] <- 1
    }
  } else {
    M[D > 0 & D <= cutoff_km] <- 1
    ## guarantee at least one neighbour (matters for Hawaii and Alaska)
    for (i in seq_len(n)) if (sum(M[i, ]) == 0) {
      o <- order(D[i, ]); M[i, o[o != i][1]] <- 1
    }
  }
  M <- pmax(M, t(M))                                   # symmetrise

  ## decay
  W <- matrix(0, n, n, dimnames = list(states, states))
  sel <- which(M > 0)
  W[sel] <- switch(decay,
    exp      = exp(-D[sel] / rho_km),
    inverse  = 1 / pmax(D[sel] / rho_km, 1e-6),
    gaussian = exp(-(D[sel] / rho_km)^2))
  W <- (W + t(W)) / 2                                  # numerical symmetry

  ## Scale so the mean row sum is 1. This does NOT change the ICAR (its
  ## precision is separately parameterised) but keeps the numbers comparable
  ## with the binary graph and helps the sampler.
  if (standardise) W <- W / mean(rowSums(W))

  if (any(rowSums(W) <= 0)) stop("isolated node(s) after applying the decay")
  if (!is_connected_graph(W))
    warning("the distance graph is DISCONNECTED; the ICAR prior assumes a ",
            "connected graph. Increase k or cutoff_km.")
  W
}


#' Is the graph connected? (breadth-first search)
is_connected_graph <- function(W) {
  n <- nrow(W); seen <- rep(FALSE, n); queue <- 1L; seen[1] <- TRUE
  while (length(queue)) {
    i <- queue[1]; queue <- queue[-1]
    nb <- which(W[i, ] > 0 & !seen)
    seen[nb] <- TRUE; queue <- c(queue, nb)
  }
  all(seen)
}


## ---------------------------------------------------------------------------
## mobility / connectivity matrix
## ---------------------------------------------------------------------------
#' Mobility matrix between states.
#'
#' @param states ordered state names
#' @param file   optional CSV with columns origin,destination,flow (state NAMES).
#'   Real sources: Census LODES / ACS commuting flows, BTS inter-state travel,
#'   SafeGraph or Meta "Data for Good" movement. Supplying real flows is
#'   strongly preferred; the gravity fallback is only a placeholder.
#' @param pop    named vector of state populations (needed for the gravity proxy)
#' @param centroids data.frame(state, lat, lon)
#' @return row-normalised matrix M (diagonal zero)
build_mobility <- function(states, file = NULL, pop = NULL, centroids = NULL,
                           fallback = c("gravity", "adjacency"),
                           adj_W = NULL) {
  fallback <- match.arg(fallback)
  n <- length(states)
  if (!is.null(file) && file.exists(file)) {
    fl <- readr::read_csv(file, show_col_types = FALSE)
    M <- matrix(0, n, n, dimnames = list(states, states))
    for (i in seq_len(nrow(fl))) {
      o <- fl$origin[i]; d <- fl$destination[i]
      if (o %in% states && d %in% states) M[o, d] <- fl$flow[i]
    }
    message("mobility: using supplied flows from ", file)
  } else if (fallback == "adjacency") {
    ## Binary land adjacency: M_ij = 1 for a shared border, 0 otherwise. The
    ## simplest possible connectivity assumption, and the easiest to describe --
    ## the covariate becomes the mean RSV rate among bordering states last week.
    ## Alaska and Hawaii share no land border, so they are joined to Washington
    ## and California, matching build_adjacency().
    message("mobility: NO flow file -> using BINARY ADJACENCY ",
            "(shared land border; AK->WA, HI->CA)")
    M <- if (!is.null(adj_W)) {
      adj_W[states, states, drop = FALSE] * 1
    } else build_adjacency(states) * 1
    M[M > 0] <- 1

  } else {
    if (is.null(pop) || is.null(centroids))
      stop("gravity proxy needs `pop` and `centroids`")
    message("mobility: NO flow file supplied -> using a GRAVITY PROXY ",
            "(pop_i * pop_j / d_ij^2). Replace with real flows when available.")
    cc <- centroids[match(states, centroids$state), ]
    R <- 6371
    rad <- function(x) x * pi / 180
    M <- matrix(0, n, n, dimnames = list(states, states))
    for (i in seq_len(n)) for (j in seq_len(n)) {
      if (i == j) next
      dlat <- rad(cc$lat[j] - cc$lat[i]); dlon <- rad(cc$lon[j] - cc$lon[i])
      a <- sin(dlat/2)^2 + cos(rad(cc$lat[i]))*cos(rad(cc$lat[j]))*sin(dlon/2)^2
      dkm <- 2 * R * asin(pmin(1, sqrt(a)))
      M[i, j] <- pop[states[i]] * pop[states[j]] / max(dkm, 50)^2
    }
  }
  diag(M) <- 0
  rs <- rowSums(M); rs[rs == 0] <- 1
  out <- M / rs                                        # row-normalised
  attr(out, "source") <- if (!is.null(file) && file.exists(file)) "flows"
                         else fallback
  attr(out, "raw") <- M                                # pre-normalisation
  out
}

## ---------------------------------------------------------------------------
## state centroids from the bundled GeoJSON
## ---------------------------------------------------------------------------
state_centroids <- function(geojson = "assets/us-states.json") {

  ## Known centroids used to FILL GAPS whatever the source. Many public
  ## us-states.json files contain only the 50 states and omit the District of
  ## Columbia; state.center likewise covers the lower 48 and places Alaska and
  ## Hawaii at fictitious positions. Any state missing after the primary source
  ## is filled from here, so the spatial graph never has a hole.
  .known <- data.frame(
    state = c("District of Columbia", "Alaska", "Hawaii"),
    lon   = c(-77.03, -152.40, -155.50),
    lat   = c( 38.90,   64.00,   19.60),
    stringsAsFactors = FALSE)

  .fill <- function(cc) {
    need <- setdiff(.known$state, cc$state)
    if (length(need)) {
      add <- .known[.known$state %in% need, ]
      message("state_centroids: filling missing centroid(s): ",
              paste(need, collapse = ", "))
      cc <- rbind(cc[, c("state", "lon", "lat")], add[, c("state", "lon", "lat")])
    }
    cc
  }

  ## --- 1. preferred: the supplied GeoJSON -----------------------------------
  if (!is.null(geojson) && length(geojson) == 1 && file.exists(geojson)) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      g <- try(jsonlite::fromJSON(geojson, simplifyVector = FALSE), silent = TRUE)
      if (!inherits(g, "try-error") && !is.null(g$features)) {
        out <- lapply(g$features, function(f) {
          co <- f$geometry$coordinates
          ring <- if (identical(f$geometry$type, "Polygon")) co[[1]] else
            co[[which.max(vapply(co, function(p) length(p[[1]]), 1L))]][[1]]
          m <- do.call(rbind, lapply(ring, function(p) c(p[[1]], p[[2]])))
          data.frame(state = f$properties$name,
                     lon = mean(m[, 1]), lat = mean(m[, 2]),
                     stringsAsFactors = FALSE)
        })
        return(.fill(do.call(rbind, out)))
      }
    }
    message("state_centroids: could not parse ", geojson, " -> built-in centres")
  } else {
    message("state_centroids: ", if (is.null(geojson)) "no geojson" else
            paste0("'", geojson, "' not found"), " -> built-in state.center")
  }

  ## --- 2. fallback: R's built-in state centres ------------------------------
  cc <- data.frame(state = datasets::state.name,
                   lon   = datasets::state.center$x,
                   lat   = datasets::state.center$y,
                   stringsAsFactors = FALSE)
  ## state.center places AK and HI at fictitious locations; drop and refill
  cc <- cc[!cc$state %in% c("Alaska", "Hawaii"), ]
  .fill(cc)
}

## ---------------------------------------------------------------------------
## main preparation
## ---------------------------------------------------------------------------
#' @param age            Epic age band ("Total", "<1 Years", ...)
#' @param start,end      analysis window
#' @param pop_file       CSV with columns state,population (e.g. from tidycensus)
#' @param climate_file   CSV with columns state,week,temp,humidity (e.g. NOAA)
#' @param mobility_file  CSV with columns origin,destination,flow
#' @param fetch_covariates if TRUE, download real population / climate / mobility
#'   via stcar_covariates.R (cached in data_cov/). Explicit *_file arguments
#'   always take precedence.
stcar_prepare <- function(age = "Total",
                          start = "2018-07-01", end = NULL,
                          pop_file = NULL, climate_file = NULL,
                          mobility_file = NULL,
                          fetch_covariates = TRUE,
                          climate_lags = 0:4,
                          mob_lags = 1:8,
                          spatial_weights = c("distance", "adjacency"),
                          mobility_fallback = c("gravity", "adjacency"),
                          knn_or_cutoff = "knn", knn_k = 6,
                          cutoff_km = 1200, decay = "exp", rho_km = 600,
                          geojson = "assets/us-states.json") {
  spatial_weights   <- match.arg(spatial_weights)
  mobility_fallback <- match.arg(mobility_fallback)

  ## centroids first: the climate fetch needs them, and this also lets the
  ## built-in fallback kick in before anything else runs
  cents <- state_centroids(geojson)

  ## ---- obtain real covariates unless files were supplied explicitly ----
  if (fetch_covariates) {
    got <- fetch_all_covariates(start = start, end = end, centroids = cents)
    if (is.null(pop_file))      pop_file      <- got$pop_file
    if (is.null(climate_file))  climate_file  <- got$climate_file
    if (is.null(mobility_file)) mobility_file <- got$mobility_file
  }

  raw <- pophive_fetch("epic")
  d <- raw |>
    mutate(week = as.Date(time),
           fips = suppressWarnings(as.integer(geography))) |>
    filter(!is.na(week), !is.na(fips), fips != 0) |>
    mutate(state = unname(FIPS_STATE[as.character(fips)])) |>
    filter(!is.na(state), age == !!age, week >= as.Date(start)) |>
    { \(x) if (is.null(end)) x else filter(x, week <= as.Date(end)) }() |>
    transmute(state, week,
              y_raw      = as.integer(epic_n_rsv),
              censored   = epic_suppressed_flag_rsv == 1,
              encounters = as.numeric(epic_n_all_encounters_weekly)) |>
    filter(encounters > 0) |>
    arrange(state, week)

  ## y is NA where censored; the interval [1,10] is passed to the likelihood
  d <- d |>
    mutate(y = ifelse(censored, NA_integer_, y_raw),
           lo = ifelse(censored, CENSOR_LO, y_raw),
           hi = ifelse(censored, CENSOR_HI, y_raw))

  states <- sort(unique(d$state))
  weeks  <- sort(unique(d$week))
  d <- d |>
    mutate(i = match(state, states),          # space index
           t = match(week, weeks))            # time index

  ## ---- population ----
  if (!is.null(pop_file) && file.exists(pop_file)) {
    pp <- readr::read_csv(pop_file, show_col_types = FALSE)
    pop <- setNames(pp$population, pp$state)[states]
  } else {
    message("!! NOTE: no population file -> using median weekly encounters as a ",
            "crude size proxy. Run fetch_population() (needs a free Census API ",
            "key) before making any population-based claim.")
    pop <- d |> group_by(state) |>
      summarise(p = median(encounters), .groups = "drop") |>
      { \(x) setNames(x$p, x$state) }()
    pop <- pop[states]
  }

  ## ---- climate ----
  ## seasonal harmonics are ALWAYS included (they are calendar, not climate)
  d <- d |> mutate(
    woy  = as.numeric(strftime(week, "%V")),
    sin1 = sin(2*pi*woy/52), cos1 = cos(2*pi*woy/52),
    sin2 = sin(4*pi*woy/52), cos2 = cos(4*pi*woy/52))

  if (!is.null(climate_file) && file.exists(climate_file)) {
    cl <- readr::read_csv(climate_file, show_col_types = FALSE) |>
      mutate(week = as.Date(week))
    d <- left_join(d, cl, by = c("state", "week"))
    ## ---- build every lag the analysis might use -------------------------
    ## Columns are named <variable>_lag<k>, including lag 0 (contemporaneous),
    ## so any subset can be selected later simply by naming it in `climate`:
    ##
    ##   CLIMATE <- c("abs_humidity_lag0", "temp_lag0")            # same week
    ##   CLIMATE <- c("abs_humidity_lag2", "temp_lag2")            # two weeks
    ##   CLIMATE <- c("abs_humidity_lag0", "abs_humidity_lag1",
    ##                "abs_humidity_lag2", "temp_lag0",
    ##                "temp_lag1", "temp_lag2")                    # several
    ##
    ## A name that was never built is silently dropped by stcar_fit_inla3(),
    ## so `climate_lags` must cover whatever the driver asks for. Lagging is
    ## done WITHIN state, after ordering by week.
    cl_vars <- intersect(c("abs_humidity", "temp", "dew_point", "precip"),
                         names(d))
    d <- d |> group_by(state) |> arrange(week)
    for (v in cl_vars) for (k in climate_lags) {
      d[[paste0(v, "_lag", k)]] <- if (k == 0) d[[v]] else dplyr::lag(d[[v]], k)
    }
    d <- ungroup(d)
    message("climate: using ", climate_file, "  (", paste(cl_vars, collapse = ", "),
            "; lags ", paste(range(climate_lags), collapse = "-"), ")")

    ## NOTE ON CHOOSING A LAG. Several lags of the same variable are highly
    ## correlated week to week, so including many of them splits one effect
    ## across several coefficients and none looks credible. Prefer choosing a
    ## lag from the empirical profile -- see lag_profile() in
    ## climate_diagnostics.R -- over entering the whole set.
  } else {
    message("!! NOTE: no climate file -> only seasonal harmonics available. ",
            "Run fetch_climate() (Open-Meteo, no API key) before making any ",
            "climatic claim.")
  }

  ## ---- graphs ----
  W <- if (identical(spatial_weights, "adjacency")) {
    build_adjacency(states)
  } else {
    build_distance_W(states, cents, type = knn_or_cutoff, k = knn_k,
                     cutoff_km = cutoff_km, decay = decay, rho_km = rho_km)
  }
  message("spatial graph: ", spatial_weights,
          " | mean neighbours = ", sprintf("%.1f", mean(rowSums(W > 0))))
  M <- build_mobility(states, file = mobility_file, pop = pop, centroids = cents,
                      fallback = mobility_fallback, adj_W = W)

  ## mobility-weighted neighbourhood RSV rate, lagged one week:
  ## this is the "imported infection pressure" covariate.
  rate <- d |>
    mutate(r = (pmax(y_raw, 1)) / encounters) |>          # crude rate for the lag
    select(state, week, r) |>
    pivot_wider(names_from = state, values_from = r) |>
    arrange(week)
  Rmat <- as.matrix(rate[, states, drop = FALSE])
  Rmat[is.na(Rmat)] <- 0
  imported <- Rmat %*% t(M)                              # (weeks x states)
  imp_df <- as.data.frame(imported); names(imp_df) <- states
  imp_df$week <- rate$week
  ## SEVERAL LAGS, not just one. A single lag of 1 week is fine for a fitting
  ## model, where every week's covariate is known. It is NOT usable for
  ## forecasting: predicting week t+h from an origin at t, a covariate lagged
  ## one week would be built from other states' counts at t+h-1, which have not
  ## happened yet. The lag must be at least as long as the horizon, so the
  ## columns are built out to `mob_lags` and the nowcast picks one that matches.
  imp_long <- imp_df |>
    pivot_longer(-week, names_to = "state", values_to = "mob_pressure") |>
    group_by(state) |>
    arrange(week)
  for (L in mob_lags)
    imp_long[[paste0("mob_pressure_lag", L)]] <- dplyr::lag(imp_long$mob_pressure, L)
  imp_long <- ungroup(imp_long)
  d <- left_join(d, imp_long, by = c("state", "week"))

  ## ---- offset ----
  d <- d |> mutate(log_offset = log(encounters),
                   log_pop    = log(unname(pop[state])))

  list(data = d, states = states, weeks = weeks,
       W = W, M = M, pop = pop, centroids = cents,
       n_state = length(states), n_week = length(weeks),
       censor = c(lo = CENSOR_LO, hi = CENSOR_HI))
}

## ---------------------------------------------------------------------------
## quick description of what was built
## ---------------------------------------------------------------------------
stcar_describe <- function(dat) {
  d <- dat$data
  cat("\n================ STCAR DATA ================\n")
  cat("states: ", dat$n_state, " | weeks: ", dat$n_week,
      " | rows: ", nrow(d), "\n", sep = "")
  cat("period: ", format(min(d$week)), " .. ", format(max(d$week)), "\n", sep = "")
  cat("censored (suppressed) rows: ",
      sprintf("%.1f%%", 100*mean(d$censored)), "\n", sep = "")
  cat("observed counts: median ", median(d$y, na.rm = TRUE),
      ", max ", max(d$y, na.rm = TRUE), "\n", sep = "")
  cat("offset log(encounters): range ",
      sprintf("%.1f", min(d$log_offset)), " .. ",
      sprintf("%.1f", max(d$log_offset)), "\n", sep = "")
  cat("adjacency: mean neighbours = ",
      sprintf("%.1f", mean(rowSums(dat$W))), "\n", sep = "")
  cat("censoring by state (worst 5):\n")
  print(d |> group_by(state) |>
          summarise(pct_censored = round(100*mean(censored), 1), .groups = "drop") |>
          arrange(desc(pct_censored)) |> head(5))
  invisible(NULL)
}


## ---------------------------------------------------------------------------
## a panel from somewhere else
## ---------------------------------------------------------------------------
#' Build the derived columns `stcar_fit_inla3()` needs, for a panel that did
#' not come from `stcar_prepare()`.
#'
#' `stcar_prepare()` does two separate jobs: it FETCHES the Epic data, and it
#' DERIVES the columns the fitters read -- the state and week indices, the
#' seasonal harmonics, the season label, the population proxy. Any other panel
#' needs the second job without the first. The CDC hub's NHSN series is the
#' case in point: it arrives already clean, but a fitter handed it directly
#' fails on every origin because `i`, `t`, `season` and `sin1` are absent.
#'
#' The failure is quiet, which is worse than its being common: `nowcast_one()`
#' returns NULL on a failed fit and `rolling_nowcast()` skips it, so sixty
#' origins can fail in a row and the only sign is a deferred warning that R
#' collapses to "50 or more warnings".
#'
#' @param d       data.frame with at least state, week, y_raw, encounters
#' @param states  the state ordering, which must match the rows and columns of
#'   `W`; taken from the data if not given
#' @param W       adjacency, for the spatial term
#' @param pop     named vector of populations, if one is available. Without it
#'   the median of `encounters` is used as a size proxy, exactly as
#'   `stcar_prepare()` does, and any population-based claim is unsupported.
stcar_panel <- function(d, states = NULL, W = NULL, centroids = NULL,
                        pop = NULL, verbose = TRUE) {
  suppressPackageStartupMessages({ library(dplyr) })

  need <- c("state", "week", "y_raw", "encounters")
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop("the panel needs: ", paste(miss, collapse = ", "), call. = FALSE)

  states <- states %||% sort(unique(d$state))
  bad <- setdiff(unique(d$state), states)
  if (length(bad))
    stop("states in the data but not in `states`: ",
         paste(head(bad, 5), collapse = ", "), call. = FALSE)

  weeks <- sort(unique(d$week))
  d <- d |>
    mutate(week = as.Date(week),
           i  = match(state, states),
           t  = match(week, weeks),
           i2 = i, t2 = t, i_m = i, t_m = t,
           woy  = as.numeric(strftime(week, "%V")),
           sin1 = sin(2 * pi * woy / 52), cos1 = cos(2 * pi * woy / 52),
           sin2 = sin(4 * pi * woy / 52), cos2 = cos(4 * pi * woy / 52),
           season = season_of(week),
           i_s1 = i, i_c1 = i, i_s2 = i, i_c2 = i)

  if (!"censored" %in% names(d)) d$censored <- FALSE
  if (!"log_offset" %in% names(d))
    d$log_offset <- log(pmax(d$encounters, 1))

  ## population: real if supplied, a size proxy otherwise. The distinction is
  ## kept visible because a proxy cannot support a per-capita claim.
  if (is.null(pop)) {
    pop <- d |> group_by(state) |>
      summarise(p = median(encounters), .groups = "drop") |>
      { \(x) setNames(x$p, x$state) }()
    if (verbose)
      message("no population supplied; using the median offset as a size ",
              "proxy. Do not interpret anything as per-capita.")
  }
  pop <- pop[states]
  d$log_pop <- log(pmax(unname(pop[d$state]), 1))

  if (is.null(W)) {
    centroids <- centroids %||% state_centroids(states)
    W <- build_distance_W(states, centroids, k = 6, decay = "exp",
                          rho_km = 600)
    if (verbose) message("built a distance graph (k = 6, rho = 600 km)")
  }
  if (!identical(dim(W), c(length(states), length(states))))
    stop("W is ", nrow(W), "x", ncol(W), " but there are ", length(states),
         " states -- the graph and the data disagree", call. = FALSE)

  if (verbose)
    cat("panel: ", nrow(d), " rows, ", length(states), " states, ",
        length(weeks), " weeks, ", sum(d$censored), " censored\n", sep = "")

  list(data = d, states = states, weeks = weeks, W = W,
       centroids = centroids %||% NULL, pop = pop)
}
