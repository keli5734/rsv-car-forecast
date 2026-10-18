## =============================================================================
## florida_seed.R — Florida as the seeding state, and distance from it
##
##   source("R/florida_seed.R")
##   dat <- add_florida_pressure(dat, weight = "exp", lag = 2)
##   res <- stcar_fit_inla3(dat, ..., extra_fixed = "fl_pressure_lag2")
##   florida_phase(res, dat)          # does the wave move outward from Florida?
##
## ---------------------------------------------------------------------------
## WHY A BARE DISTANCE VECTOR CANNOT WORK
## ---------------------------------------------------------------------------
## The suggestion is to add distance from Florida as a 51-vector covariate. As a
## MAIN EFFECT that cannot be identified: distance is constant within a state,
## and the model already has b_i, a free effect per state. Any distance
## coefficient could be absorbed into b_i without changing the fit at all, so
## the two are perfectly collinear and INLA would return an arbitrary split
## between them.
##
## The idea is still right -- it just has to enter somewhere it varies in time.
## Two forms do that, and they test different claims:
##
##  1. SEEDING PRESSURE (fitted here). Florida's CURRENT level, discounted by
##     distance:
##
##         x_it = w(d_i) * r_FL,t-L
##
##     where w falls with distance and L is a lag. This varies over t, so it is
##     identified, and it is directly usable in a forecast: at origin t it is
##     known for t+h whenever L >= h. It asks "does RSV in Florida now predict
##     RSV here later, more strongly for nearby states?"
##
##  2. PHASE GRADIENT (diagnosed, not fitted). If the epidemic radiates from
##     Florida, states further away should PEAK LATER -- their seasonal phase
##     psi_i should increase with distance. `florida_phase()` extracts psi_i
##     from a fitted harmonic model and regresses it on distance. This is the
##     cleaner test of the supervisor's rationale, because it looks at the
##     travelling wave itself rather than at a week-to-week association.
##
## Reporting both is stronger than either: the phase gradient says whether the
## wave exists, the pressure covariate says whether it is useful for prediction.
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Distance from a seeding state to every other, in km.
seed_distance <- function(states, centroids = NULL, seed = "Florida") {
  if (is.null(centroids)) {
    if (!exists("state_centroids"))
      stop("supply `centroids`, or source stcar_prep.R for state_centroids()",
           call. = FALSE)
    centroids <- state_centroids(states)
  }
  ## MATCH BY NAME, never by position. state_centroids() may return its rows in
  ## its own order, and relying on the order silently attaches every distance to
  ## the wrong state -- which is hard to notice because the numbers still look
  ## like distances.
  nmcol <- intersect(c("state", "NAME", "name"), names(centroids))[1]
  if (is.na(nmcol))
    stop("`centroids` needs a state-name column so rows can be matched by name",
         call. = FALSE)
  idx <- match(states, centroids[[nmcol]])
  if (anyNA(idx))
    stop("no centroid for: ", paste(states[is.na(idx)], collapse = ", "),
         call. = FALSE)
  centroids <- centroids[idx, , drop = FALSE]

  if (!seed %in% states) stop(seed, " is not in `states`", call. = FALSE)
  i <- match(seed, states)
  R <- 6371
  rad <- function(x) x * pi / 180
  la <- rad(centroids$lat); lo <- rad(centroids$lon)
  dlat <- la - la[i]; dlon <- lo - lo[i]
  a <- sin(dlat / 2)^2 + cos(la[i]) * cos(la) * sin(dlon / 2)^2
  d <- 2 * R * asin(pmin(1, sqrt(a)))
  out <- setNames(d, states)

  ## a cheap check that the alignment is right: the seed must be 0 from itself,
  ## and its neighbours must be closer than the far coast
  if (abs(out[[seed]]) > 1e-6)
    warning("distance from ", seed, " to itself is ", round(out[[seed]]),
            " km -- the centroids are not aligned with `states`", call. = FALSE)
  out
}


#' Add a seeding-pressure covariate: the seed state's own rate, discounted by
#' distance (or by commuting share) and lagged.
#'
#' @param weight  How the seed state's signal is discounted. These are
#'   different hypotheses about HOW the epidemic spreads, and the phase analysis
#'   suggested they are not equivalent: once latitude was accounted for,
#'   distance from Florida was no longer independently predictive of peak
#'   timing. Test them separately rather than assuming distance.
#'
#'                "exp"      w = exp(-d / rho), great-circle distance. Radial
#'                           spread outward from the seed.
#'                "inverse"  w = 1 / (1 + d / rho), a heavier tail than "exp".
#'                "rank"     w = 1 - rank(d) / n; only the ordering matters, so
#'                           the result does not depend on the decay scale.
#'                "latitude" w = exp(-(lat_i - lat_seed) / rho_deg), so the
#'                           discount depends on how far NORTH a state is rather
#'                           than how far away. This is the south-to-north wave
#'                           the RSV literature describes, and the axis the
#'                           phase regression favoured.
#'                "lat_rank" the same, by latitude ordering only.
#'                "commute"  w = the commuting share with the seed state.
#'                           Distance is a proxy for contact; commuting is a
#'                           measurement of it.
#' @param rho_deg decay scale for the latitude weights, in degrees
#' @param lag     weeks. For FORECASTING this must be at least the horizon, or
#'                the covariate on a future row uses a Florida value that has
#'                not been observed at the origin.
add_florida_pressure <- function(dat, seed = "Florida", lag = 2,
                                 weight = c("exp", "inverse", "rank",
                                            "latitude", "lat_rank", "commute"),
                                 rho_km = 1500, rho_deg = 10,
                                 suffix = NULL, exclude_seed = TRUE,
                                 verbose = TRUE) {
  weight <- match.arg(weight)
  d <- dat$data
  states <- dat$states

  dist <- seed_distance(states, dat$centroids %||% NULL, seed)

  ## latitude, matched by name for the same reason distances are
  cen <- dat$centroids %||% state_centroids(states)
  nmcol <- intersect(c("state", "NAME", "name"), names(cen))[1]
  lat <- setNames(cen$lat[match(states, cen[[nmcol]])], states)
  dlat <- lat - lat[[seed]]        # degrees north of the seed

  w <- switch(weight,
    exp      = exp(-dist / rho_km),
    inverse  = 1 / (1 + dist / rho_km),
    rank     = 1 - (rank(dist) - 1) / (length(dist) - 1),
    latitude = exp(-pmax(dlat, 0) / rho_deg),
    lat_rank = 1 - (rank(dlat) - 1) / (length(dlat) - 1),
    commute = {
      M <- dat$M
      if (is.null(M) || !seed %in% rownames(M))
        stop("no commuting matrix with ", seed, " in dat$M; ",
             "rebuild with stcar_prepare(fetch_covariates = TRUE)", call. = FALSE)
      v <- M[, seed]; v / max(v)
    })
  w[seed] <- if (exclude_seed) 0 else 1      # Florida cannot seed itself

  ## the seed state's own rate, by week
  fl <- d |> filter(state == seed) |>
    transmute(week, r_seed = ifelse(censored, NA_real_,
                                    y_raw / pmax(encounters, 1))) |>
    arrange(week)
  ## a suppressed seed week leaves a gap; carrying the last value forward is
  ## the least-bad fill, and those weeks are the off-season where the
  ## covariate matters least
  fl$r_seed <- zoo::na.locf(fl$r_seed, na.rm = FALSE)
  fl$r_seed <- zoo::na.locf(fl$r_seed, fromLast = TRUE, na.rm = FALSE)

  ## the weight is in the column name, so two weightings can sit side by side
  ## in the same data and be compared without refitting the other
  col <- paste0("fl_", suffix %||% weight, "_lag", lag)
  fl[[col]] <- dplyr::lag(fl$r_seed, lag)
  fl[[col]] <- as.numeric(scale(fl[[col]]))

  ## Drop any previous version of this column first. Without it, calling the
  ## function twice makes the join produce `<col>.x` and `<col>.y` and the bare
  ## name disappears -- which surfaces as a confusing "column not found" inside
  ## the mutate rather than at the join.
  if (col %in% names(d)) {
    message("replacing the existing ", col, " column")
    d[[col]] <- NULL
  }

  d <- d |>
    left_join(fl |> select(week, all_of(col)), by = "week") |>
    mutate(!!col := .data[[col]] * unname(w[state]))

  if (verbose) {
    cat("\n=========== SEEDING PRESSURE FROM ", toupper(seed),
        " ===========\n", sep = "")
    cat("weight: ", weight,
        if (weight %in% c("exp", "inverse")) paste0(" (rho = ", rho_km, " km)") else "",
        ";  lag: ", lag, " weeks\n", sep = "")
    cat("column: ", col, "\n", sep = "")
    ex <- data.frame(state = states, km = round(dist),
                     deg_north = round(dlat, 1), w = round(w, 3))
    ex <- ex[order(-ex$w), ]
    cat("\nhighest and lowest weights:\n")
    print(rbind(head(ex[ex$state != seed, ], 5), tail(ex, 3)), row.names = FALSE)
    cat("\nNOTE: distance alone is collinear with the state effect b_i and\n",
        "cannot be fitted as a main effect. This column is distance TIMES the\n",
        "seed state's lagged rate, so it varies in time and is identified.\n",
        sep = "")
  }
  dat$data <- d
  attr(dat, "seed_weight") <- w
  attr(dat, "seed_distance") <- dist
  dat
}


## ---------------------------------------------------------------------------
## the phase gradient
## ---------------------------------------------------------------------------
#' Do states further from Florida peak later?
#'
#' For a harmonic model the first pair of coefficients gives each state an
#' amplitude and a phase:
#'
#'     a_i sin + c_i cos = A_i sin(2 pi w / 52 + psi_i)
#'
#' and psi_i converted to weeks is that state's peak timing. If the epidemic
#' radiates from Florida, peak week should rise with distance. This is the
#' direct test of the seeding hypothesis and does not require fitting any
#' covariate at all -- it reads the answer out of a model already fitted.
florida_phase <- function(res, dat, seed = "Florida",
                          out_dir = "out/florida", base = 13) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rnd <- res$fit$summary.random
  nm <- names(rnd)
  si <- grep("i_s1$", nm, value = TRUE)[1]
  ci <- grep("i_c1$", nm, value = TRUE)[1]
  if (is.na(si) || is.na(ci))
    stop("no harmonic random slopes in this fit; refit with ",
         "interaction = 'harmonic', 'harmonic_ar1' or 'harmonic_season'",
         call. = FALSE)

  a <- rnd[[si]]$mean; cc <- rnd[[ci]]$mean
  ## the random-effect rows are indexed by dat$states order, and the ID column
  ## is carried so the mapping can be checked rather than assumed
  ids <- rnd[[si]]$ID
  states <- dat$states[ids]
  amp <- sqrt(a^2 + cc^2)

  ## PEAK TIMING. The harmonics are built on WEEK OF YEAR, so a peak in late
  ## December lands near 52 and one in early January near 1 -- the same
  ## calendar moment at opposite ends of the axis. Regressing that directly is
  ## meaningless, and produces exactly the null-with-no-fit seen before.
  ##
  ## Converting to WEEK OF SEASON fixes it: with the season starting in July,
  ## a December/January peak sits near week 26, in the middle of the range,
  ## and no state wraps.
  psi <- atan2(cc, a)
  peak_woy <- ((pi / 2 - psi) / (2 * pi) * 52) %% 52
  peak_week <- (peak_woy - 26) %% 52          # July = week 0

  dist <- seed_distance(dat$states, dat$centroids %||% NULL, seed)
  tab <- data.frame(state = states, km = unname(dist[states]),
                    amplitude = amp,
                    peak_woy = peak_woy, peak_week = peak_week) |>
    arrange(km)

  cat("\npeak timing (week of season, July = 0):  median ",
      sprintf("%.1f", median(tab$peak_week)), ", range ",
      sprintf("%.1f", min(tab$peak_week)), " to ",
      sprintf("%.1f", max(tab$peak_week)), "\n", sep = "")
  if (diff(range(tab$peak_week)) > 40)
    cat("  !! the range still spans most of a year, so some states are wrapping.\n",
        "  Peak timing is not comparable across them and the regression below\n",
        "  should not be read.\n", sep = "")

  ## These harmonics are averaged over EVERY season in the fit, including the
  ## pandemic years when RSV was displaced by months. A state whose 2021 season
  ## arrived in summer has its average phase pulled away from its usual one, so
  ## a genuine gradient would be blurred. Refitting on 2022 onward, or using
  ## interaction = "harmonic_season" and reading phase per season, is the
  ## sharper test.
  fit <- lm(peak_week ~ km, data = tab)
  s <- summary(fit)
  b <- coef(s)["km", ]
  cat("\n=========== PHASE GRADIENT FROM ", toupper(seed), " ===========\n", sep = "")
  cat("peak week vs distance:  slope = ", sprintf("%.3f", b[1] * 1000),
      " weeks per 1000 km  (SE ", sprintf("%.3f", b[2] * 1000),
      ", p = ", format.pval(b[4], digits = 2), ")\n", sep = "")
  cat("R-squared: ", sprintf("%.3f", s$r.squared), "\n", sep = "")
  cat(if (b[4] < 0.05 && b[1] > 0)
    "  peaks are later further from the seed state, consistent with a wave\n  moving outward.\n"
    else if (b[4] < 0.05 && b[1] < 0)
      "  peaks are EARLIER further away, which is the opposite of the seeding\n  hypothesis and needs explaining.\n"
    else "  no credible gradient: these data do not show peak timing varying\n  with distance from the seed state.\n")
  write_csv(tab, file.path(out_dir, "phase_by_distance.csv"))

  p <- ggplot(tab, aes(km, peak_week)) +
    geom_smooth(method = "lm", se = TRUE, colour = "#b2182b", fill = "#b2182b",
                alpha = .12, linewidth = .7, formula = y ~ x) +
    geom_point(aes(size = amplitude), alpha = .7, colour = "#2166ac") +
    ggrepel::geom_text_repel(aes(label = state), size = 2.4, max.overlaps = 12,
                             colour = "grey35") +
    scale_size_continuous(range = c(1.5, 5), name = "seasonal\namplitude") +
    labs(title = paste0("Does the season start in ", seed, " and move outward?"),
         subtitle = sprintf("slope %.2f weeks per 1000 km (p = %s)",
                            b[1] * 1000, format.pval(b[4], digits = 2)),
         x = paste0("distance from ", seed, " (km)"),
         y = "fitted peak week of season (July = 0)") +
    theme_minimal(base_size = base)
  ggsave(file.path(out_dir, "phase_by_distance.png"), p,
         width = 8, height = 5.5, dpi = 300)
  invisible(list(table = tab, fit = fit, plot = p))
}


## ---------------------------------------------------------------------------
## phase season by season
## ---------------------------------------------------------------------------
#' Peak timing per state AND season, from a harmonic_season fit.
#'
#' `florida_phase()` reads one phase per state, averaged over every season in
#' the fit. That answers "do states typically peak later further north", but it
#' cannot distinguish two very different situations:
#'
#'   a TRAVELLING WAVE, where the same ordering of states recurs every year --
#'   strong evidence that the epidemic propagates geographically; and
#'
#'   a STABLE AVERAGE over seasons whose orderings differ, which would produce
#'   the same mean gradient while the wave itself is not reproducible.
#'
#' Fitting the gradient WITHIN each season separates them. It also sidesteps the
#' pandemic problem: rather than averaging 2020-21 and 2021-22 into every
#' state's phase, those seasons get their own rows and can be reported apart.
#'
#' Requires interaction = "harmonic_season", which carries both a state-level
#' harmonic pair and a state-season one. The state-season deviation is added to
#' the state mean, so psi_is is that season's own phase for that state.
florida_phase_by_season <- function(res, dat, seed = "Florida",
                                    min_amp = 0.1, min_states = 20,
                                    out_dir = "out/florida", base = 13) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rnd <- res$fit$summary.random; nm <- names(rnd)
  s_st <- grep("^i_s1$", nm, value = TRUE)[1]
  c_st <- grep("^i_c1$", nm, value = TRUE)[1]
  s_is <- grep("^is_s1$", nm, value = TRUE)[1]
  c_is <- grep("^is_c1$", nm, value = TRUE)[1]
  if (any(is.na(c(s_st, c_st, s_is, c_is))))
    stop("this fit has no state-season harmonics; refit with ",
         "interaction = 'harmonic_season'", call. = FALSE)

  ## the state-season index was built as factor(paste(state, season)); rebuild
  ## the same mapping so the coefficients can be labelled
  d <- res$data
  if (!"season" %in% names(d)) d$season <- season_of(d$week)
  key <- sort(unique(paste(d$state, d$season)))
  lab <- data.frame(idx = seq_along(key),
                    state = sub(" [0-9]{4}-[0-9]{4}$", "", key),
                    season = sub("^.* ", "", key), stringsAsFactors = FALSE)

  st <- data.frame(state = dat$states[rnd[[s_st]]$ID],
                   a_st = rnd[[s_st]]$mean, c_st = rnd[[c_st]]$mean)
  is <- data.frame(idx = rnd[[s_is]]$ID,
                   a_is = rnd[[s_is]]$mean, c_is = rnd[[c_is]]$mean) |>
    left_join(lab, by = "idx")

  ## the season's own coefficients are the state mean plus its deviation
  tab <- is |>
    left_join(st, by = "state") |>
    mutate(a = a_st + a_is, cc = c_st + c_is,
           amplitude = sqrt(a^2 + cc^2),
           peak_woy = ((pi / 2 - atan2(cc, a)) / (2 * pi) * 52) %% 52,
           peak_week = (peak_woy - 26) %% 52)          # July = 0

  dist <- seed_distance(dat$states, dat$centroids %||% NULL, seed)
  cen <- dat$centroids %||% state_centroids(dat$states)
  nmc <- intersect(c("state", "NAME", "name"), names(cen))[1]
  tab$km  <- unname(dist[tab$state])
  tab$lat <- cen$lat[match(tab$state, cen[[nmc]])]

  n_flat <- sum(tab$amplitude < min_amp)
  if (n_flat)
    message("dropping ", n_flat, " state-seasons with no identifiable seasonal ",
            "signal (amplitude < ", min_amp, ")")
  tab <- filter(tab, amplitude >= min_amp)

  ## ---- the gradient within each season ----------------------------------
  ## Seasons with few usable states cannot support a regression; they are
  ## reported with their n rather than dropped silently, because a season being
  ## unusable is itself informative about how much the data supports.
  fits <- tab |> group_by(season) |> group_modify(function(g, k) {
    if (nrow(g) < min_states || diff(range(g$peak_week)) > 40)
      return(data.frame(n = nrow(g), slope_km = NA_real_, p_km = NA_real_,
                        slope_lat = NA_real_, p_lat = NA_real_,
                        R2 = NA_real_, note = if (nrow(g) < min_states)
                          "too few states" else "phases wrap"))
    m1 <- summary(lm(peak_week ~ km, data = g))
    m2 <- summary(lm(peak_week ~ lat, data = g))
    data.frame(n = nrow(g),
               slope_km = m1$coefficients["km", 1] * 1000,
               p_km = m1$coefficients["km", 4],
               slope_lat = m2$coefficients["lat", 1],
               p_lat = m2$coefficients["lat", 4],
               R2 = m1$r.squared, note = "")
  }) |> ungroup()

  cat("\n=========== PHASE GRADIENT, SEASON BY SEASON ===========\n")
  cat("slope_km is weeks per 1000 km; slope_lat is weeks per degree north\n\n")
  print(as.data.frame(fits), row.names = FALSE, digits = 3)

  ok <- fits |> filter(!is.na(slope_lat))
  if (nrow(ok)) {
    cat("\nseasons with a usable fit: ", nrow(ok), "\n", sep = "")
    cat("latitude slope positive in ", sum(ok$slope_lat > 0), " of them",
        ", credibly so (p < 0.05) in ", sum(ok$p_lat < 0.05), "\n", sep = "")
    cat(if (sum(ok$slope_lat > 0) == nrow(ok))
      "  the gradient points the same way in EVERY season, which is what a\n  travelling wave looks like.\n"
      else if (sum(ok$slope_lat > 0) >= 0.75 * nrow(ok))
        "  the gradient points the same way in most seasons; the exceptions are\n  worth naming rather than averaged away.\n"
      else "  the direction is not consistent across seasons. An average gradient\n  can then arise without a reproducible wave, and should not be\n  described as one.\n")
  }

  ## ---- is the ORDERING of states stable? --------------------------------
  ## The sharper test. A wave means the same states lead and the same states
  ## follow, year after year. Correlating each season's ranking against the
  ## others measures that directly, and does not depend on distance or latitude
  ## being the right axis.
  wide <- tab |> select(state, season, peak_week) |>
    tidyr::pivot_wider(names_from = season, values_from = peak_week)
  mm <- as.matrix(wide[, -1, drop = FALSE])
  keep <- colSums(!is.na(mm)) >= min_states
  mm <- mm[, keep, drop = FALSE]
  if (ncol(mm) >= 2) {
    cr <- suppressWarnings(cor(mm, use = "pairwise.complete.obs",
                               method = "spearman"))
    off <- cr[upper.tri(cr)]
    cat("\nagreement between seasons in WHICH states peak first\n")
    cat("  median rank correlation: ", sprintf("%.2f", median(off, na.rm = TRUE)),
        ",  range ", sprintf("%.2f", min(off, na.rm = TRUE)), " to ",
        sprintf("%.2f", max(off, na.rm = TRUE)), "\n", sep = "")
    cat(if (median(off, na.rm = TRUE) > 0.5)
      "  the ordering largely recurs, which is stronger evidence for a wave\n  than the average gradient alone.\n"
      else "  the ordering does not recur. Peak timing differs by state on\n  average without the same states leading each year.\n")
    write.csv(cr, file.path(out_dir, "phase_season_agreement.csv"))
  }

  readr::write_csv(tab, file.path(out_dir, "phase_by_state_season.csv"))
  readr::write_csv(fits, file.path(out_dir, "phase_gradient_by_season.csv"))

  p <- ggplot(tab, aes(lat, peak_week)) +
    geom_smooth(method = "lm", se = FALSE, colour = "#b2182b",
                linewidth = .6, formula = y ~ x) +
    geom_point(alpha = .6, size = 1.6, colour = "#2166ac") +
    facet_wrap(~ season) +
    labs(title = "Peak timing against latitude, season by season",
         subtitle = paste("each point is one state;  a wave would tilt the same",
                          "way every year"),
         x = "latitude", y = "peak week of season (July = 0)") +
    theme_minimal(base_size = base) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(out_dir, "phase_by_season.png"), p,
         width = 10, height = 7, dpi = 300)

  invisible(list(table = tab, by_season = fits, plot = p))
}
