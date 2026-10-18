 
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readr); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a


## ---------------------------------------------------------------------------
## distributional scores, used by both scorers
## ---------------------------------------------------------------------------
## Defined once at file scope. Duplicating them inside each function invites the
## two copies to drift apart, and a scoring rule that differs between the
## pooled and per-state tables would be very hard to notice.

#' Weighted interval score with a single central interval.
#'
#' CRPS cannot be computed from a median and one interval alone; WIS is the
#' quantile-based approximation used by the COVID and influenza Forecast Hubs:
#'
#'   IS  = (u - l) + (2/alpha)(l - y)1{y<l} + (2/alpha)(y - u)1{y>u}
#'   WIS = [ |y - m|/2 + (alpha/2) IS ] / 1.5
#'
#' It penalises a wide interval and punishes one that misses, so accuracy cannot
#' be bought with vagueness. On the scale of the counts, like MAE. Lower better.
.wis1 <- function(y, m, l, u, alpha = 0.05) {
  is <- (u - l) + (2 / alpha) * pmax(l - y, 0) + (2 / alpha) * pmax(y - u, 0)
  (0.5 * abs(y - m) + (alpha / 2) * is) / 1.5
}

#' Exact CRPS for a Poisson predictive.
#'
#'   CRPS = y(2F(y) - 1) - mu(2F(y-1) - 1) - mu e^{-2mu}(I0(2mu) + I1(2mu))
#'
#' This is a true CRPS, but for a different object: it treats the fitted mean as
#' known and so ignores the parameter uncertainty the interval represents. A
#' cross-check on the point predictive, not a substitute for WIS.
.crps_pois <- function(y, mu) {
  mu <- pmax(mu, 1e-10)
  Fy  <- stats::ppois(y, mu)
  Fy1 <- stats::ppois(y - 1, mu)
  b <- besselI(2 * mu, 0, expon.scaled = TRUE) +
       besselI(2 * mu, 1, expon.scaled = TRUE)
  y * (2 * Fy - 1) - mu * (2 * Fy1 - 1) - mu * b
}

## ---------------------------------------------------------------------------
## one origin
## ---------------------------------------------------------------------------
#' Fit at a single origin and return the predictions for t+1 .. t+max(horizons).
#'
#' @param dat        stcar_prepare() output
#' @param origin     last week the model is allowed to see
#' @param horizons   weeks ahead to predict
#' @param ...        model settings, passed to stcar_fit_inla3()
#' @param censored  whether the OBSERVED weeks carry an interval-censored
#'   likelihood. TRUE for the Epic panel, where counts at or below the ceiling
#'   are suppressed; FALSE for a fully reported series such as the CDC hub's
#'   NHSN target. The weeks after the origin are always marked uncensored and
#'   left missing regardless -- an unobserved week is not a suppressed one, and
#'   conflating the two would have the model treat every forecast week as
#'   known to lie in [0, ceiling].
nowcast_one <- function(dat, origin, horizons = 1:4,
                        climate = NULL, use_pop = FALSE, mobility = "none",
                        mob_lag = NULL, extra_fixed = NULL, censored = TRUE,
                        hyper_priors = NULL, quantiles = NULL,
                        interaction = "type2", add_trend = TRUE,
                        slope_model = "none", family = "cenpoisson",
                        cen_lo = 1, cen_hi = 10,
                        subset_states = NULL, verbose = FALSE) {

  origin <- as.Date(origin)
  d0 <- dat$data
  if (!is.null(subset_states)) d0 <- dplyr::filter(d0, state %in% subset_states)

  wks <- sort(unique(d0$week))
  targets <- wks[wks > origin][seq_len(max(horizons))]
  targets <- targets[!is.na(targets)]
  if (!length(targets)) return(NULL)

  ## hide everything after the origin
  is_fc <- d0$week %in% targets
  d0$.truth    <- ifelse(is_fc, as.numeric(d0$y_raw), NA_real_)
  d0$.forecast <- is_fc
  d0$y_raw[is_fc]    <- NA_integer_
  d0$censored[is_fc] <- FALSE          # unobserved is not suppressed

  ## weeks strictly after the last target are dropped: they carry no
  ## information for this origin and only slow the fit
  d0 <- d0 |> filter(week <= max(targets))

  dat2 <- dat; dat2$data <- d0
  ## The other-state covariate must be lagged at least as far as the longest
  ## horizon. Otherwise the value on a future row is built from counts that had
  ## not occurred at the origin, the forecast sees the future, and the scores
  ## are meaningless. Defaulting to max(horizons) makes that impossible to get
  ## wrong by accident.
  if (mobility != "none") {
    if (is.null(mob_lag)) mob_lag <- max(horizons)
    if (mob_lag < max(horizons))
      stop("mob_lag = ", mob_lag, " is shorter than the longest horizon (",
           max(horizons), "). The covariate would be built from data that is ",
           "not available at the origin.", call. = FALSE)
  }

  args <- list(dat2, climate = climate, use_pop = use_pop, mobility = mobility,
               slope_model = slope_model, family = family,
               interaction = interaction, censored = censored,
               subset_states = subset_states,
               cen_lo = cen_lo, cen_hi = cen_hi, verbose = verbose)
  if ("add_trend" %in% names(formals(stcar_fit_inla3))) args$add_trend <- add_trend
  ## the nowcast reads only fitted values, so ask for the stripped-down fit:
  ## config/dic/waic make the object several times larger for no benefit here
  if ("light" %in% names(formals(stcar_fit_inla3))) args$light <- TRUE
  ## priors borrowed from another fit, passed straight through. They are a
  ## property of the MODEL, not of the origin, so the same set is used at every
  ## origin -- nothing from the series being forecast enters them.
  if (!is.null(hyper_priors)) args$hyper_priors <- hyper_priors
  if (!is.null(quantiles)) args$quantiles <- quantiles
  if (mobility != "none") args$mob_lag <- mob_lag

  ## Any covariate built outside the fitter -- seeding pressure, say. Its LAG
  ## must be at least the longest horizon: a column named *_lagL carries the
  ## source state's value from L weeks before, so on a row at t+h it uses data
  ## from t+h-L. That is only available at the origin when L >= h. A shorter lag
  ## lets the forecast see the future and the scores become meaningless, so the
  ## lag is read out of the column name and checked.
  if (length(extra_fixed)) {
    lags <- suppressWarnings(as.integer(sub(".*_lag(\\d+)$", "\\1", extra_fixed)))
    bad <- extra_fixed[!is.na(lags) & lags < max(horizons)]
    if (length(bad))
      stop("these covariates are lagged less than the longest horizon (",
           max(horizons), "): ", paste(bad, collapse = ", "),
           "\n  they would use values not yet observed at the origin. Rebuild ",
           "them with lag = ", max(horizons), " or more.", call. = FALSE)
    if (any(is.na(lags)))
      warning("cannot read a lag from: ",
              paste(extra_fixed[is.na(lags)], collapse = ", "),
              " -- check it is knowable at the origin", call. = FALSE)
    args$extra_fixed <- extra_fixed
  }

  res <- try(do.call(stcar_fit_inla3, args), silent = TRUE)
  if (inherits(res, "try-error")) {
    warning("origin ", format(origin), " failed: ",
            conditionMessage(attr(res, "condition")), call. = FALSE)
    return(NULL)
  }

  out <- res$data
  k <- function(x) paste(x$state, x$week)
  out$.truth <- d0$.truth[match(k(out), k(d0))]
  ## keep every q* column the fit produced, so a hub submission has its 23
  ## levels without a second pass
  out <- out |>
    filter(week %in% targets) |>
    mutate(origin  = origin,
           horizon = as.numeric(difftime(week, origin, units = "weeks"))) |>
    select(state, week, origin, horizon, .truth,
           mu_med, mu_lo, mu_hi, encounters)
  rm(res); invisible(gc())
  out
}


## ---------------------------------------------------------------------------
## roll the origin through the validation season
## ---------------------------------------------------------------------------
#' @param valid_start  first origin; everything before it is training only
#' @param refit_every  refit every k weeks. 1 is exact and matches real
#'                     deployment, where a model is refitted as data arrives;
#'                     larger values save time and MUST be reported.
#' @param checkpoint_dir  each origin is written separately, so an interrupted
#'                     run resumes instead of restarting
rolling_nowcast <- function(dat, valid_start, valid_end = NULL,
                            horizons = 1:4, refit_every = 1,
                            checkpoint_dir = "out/nowcast_inla/origins",
                            resume = TRUE, accumulate = FALSE,
                            verbose = FALSE, ...) {

  wks <- sort(unique(dat$data$week))
  vs <- as.Date(valid_start)
  ve <- if (is.null(valid_end)) max(wks) - 7 * max(horizons) else as.Date(valid_end)
  origins <- wks[wks >= vs & wks <= ve]
  origins <- origins[seq(1, length(origins), by = refit_every)]

  if (!is.null(checkpoint_dir))
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  ck <- function(o) file.path(checkpoint_dir, paste0("origin_", format(o), ".rds"))

  message("rolling nowcast: ", length(origins), " origins, ",
          format(min(origins)), " to ", format(max(origins)),
          "  (refit every ", refit_every, " weeks)")

  ## With `accumulate = FALSE` (the default) nothing is kept in memory during
  ## the loop: each origin is written to its own file and they are read back at
  ## the end. Holding results for fifty origins is not itself large, but any
  ## reference kept alive can prevent gc() from releasing the fit that produced
  ## it, which is the usual reason a long run grows until it crashes.
  res <- list(); t0 <- Sys.time()
  for (i in seq_along(origins)) {
    o <- origins[i]
    if (!is.null(checkpoint_dir) && resume && file.exists(ck(o))) {
      if (accumulate) res[[length(res) + 1]] <- readRDS(ck(o))
      message("[", i, "/", length(origins), "] ", format(o), " cached"); next
    }
    message("[", i, "/", length(origins), "] ", format(o))
    r <- nowcast_one(dat, o, horizons = horizons, verbose = verbose, ...)
    ## A failed fit returns NULL. Skipping quietly is right for an occasional
    ## non-convergence, but if the FIRST origin fails the cause is almost
    ## always structural -- a missing column, a graph of the wrong size -- and
    ## every later origin will fail the same way. R defers warnings and
    ## collapses them to "50 or more warnings", so a silent skip there means
    ## sixty failures and no visible reason. The first one therefore stops.
    if (is.null(r)) {
      if (i == 1)
        stop("the first origin produced no fit, so the rest will fail the ",
             "same way.\n  The warning above gives the reason. A panel not ",
             "built by stcar_prepare()\n  usually needs stcar_panel() first: ",
             "the fitter reads i, t, season, sin1,\n  cos1 and log_pop, none ",
             "of which a raw data frame carries.", call. = FALSE)
      next
    }
    if (!is.null(checkpoint_dir)) saveRDS(r, ck(o))
    if (accumulate) res[[length(res) + 1]] <- r

    rm(r); invisible(gc())
    if (verbose || i %% 5 == 0)
      message(sprintf("    memory in use: %.2f GB",
                      sum(gc()[, "used"] * 8) / 1e6))

    if (i == 1) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message(sprintf("  first origin took %.1f min -> about %.1f h for %d origins",
                      el, el * length(origins) / 60, length(origins)))
    }
  }

  ## assemble from disk
  if (!accumulate && !is.null(checkpoint_dir)) {
    fs <- list.files(checkpoint_dir, pattern = "^origin_.*\\.rds$",
                     full.names = TRUE)
    if (!length(fs)) stop("no origin files in ", checkpoint_dir, call. = FALSE)
    message("assembling ", length(fs), " origins from disk")
    return(bind_rows(lapply(fs, readRDS)))
  }
  if (!length(res)) stop("no origins produced predictions", call. = FALSE)
  bind_rows(res)
}


## ---------------------------------------------------------------------------
## scoring
## ---------------------------------------------------------------------------
#' Score against two baselines.
#'
#' PERSISTENCE carries the last observed count forward and is hard to beat at
#' short horizons. SEASONAL NAIVE uses the same week a year earlier and is the
#' meaningful comparator for a seasonal disease.
#'
#' Suppressed target weeks have no point truth, so they are scored separately:
#' the only check available there is whether the predictive interval overlaps
#' the censoring interval.
score_nowcast_inla <- function(nc, dat, cen_hi = 10) {
  src <- dat$data |> transmute(state, week, y = as.numeric(y_raw),
                               censored)
  d <- nc |> left_join(src |> transmute(state, week, cen_t = censored),
                       by = c("state", "week"))

  ## persistence: the count at the origin
  pers <- src |> transmute(state, origin = week, persist = ifelse(censored, 5.5, y))
  ## seasonal naive: the same week one year earlier
  sn <- src |> transmute(state, w52 = week, seas = ifelse(censored, 5.5, y))

  d <- d |>
    left_join(pers, by = c("state", "origin")) |>
    mutate(w52 = week - 7 * 52) |>
    left_join(sn, by = c("state", "w52"))

  ## the distributional scores need to be computed before they are summarised
  d <- d |>
    mutate(WIS = .wis1(.truth, mu_med, mu_lo, mu_hi, 0.05),
           CRPS_pois = .crps_pois(.truth, mu_med))

  sc <- d |>
    filter(!cen_t, is.finite(.truth)) |>
    group_by(horizon) |>
    summarise(
      n        = n(),
      MAE      = mean(abs(mu_med - .truth)),
      RMSE     = sqrt(mean((mu_med - .truth)^2)),
      WIS      = mean(WIS, na.rm = TRUE),
      CRPS_pois = mean(CRPS_pois, na.rm = TRUE),
      bias     = mean(mu_med - .truth),
      corr     = suppressWarnings(cor(mu_med, .truth)),
      coverage = mean(.truth >= mu_lo & .truth <= mu_hi),
      MAE_persistence = mean(abs(persist - .truth), na.rm = TRUE),
      MAE_seasonal    = mean(abs(seas - .truth), na.rm = TRUE),
      .groups = "drop") |>
    mutate(skill_vs_persistence = 1 - MAE / MAE_persistence,
           skill_vs_seasonal    = 1 - MAE / MAE_seasonal,
           ## skill on the proper score, which is the one to lead with: it
           ## credits a sharp well-placed interval and penalises a vague one,
           ## where MAE is blind to both
           wis_skill_persistence = 1 - WIS / MAE_persistence)

  cen <- d |> filter(cen_t) |>
    group_by(horizon) |>
    summarise(n_censored = n(),
              interval_overlap = mean(mu_hi >= 1 & mu_lo <= cen_hi),
              .groups = "drop")

  out <- left_join(sc, cen, by = "horizon")
  cat("\n=========== NOWCAST SCORES ===========\n")
  print(as.data.frame(out), row.names = FALSE, digits = 3)
  ## Is the model just reproducing persistence? A random walk extrapolates
  ## flat, so its h-step forecast is close to the value at the origin, and the
  ## prediction plotted at t+h looks like the data shifted right by h weeks.
  ## A correlation with persistence above about 0.95 means the two are the same
  ## forecast and no amount of fitting will separate them.
  lag_chk <- d |> filter(!cen_t, is.finite(.truth), is.finite(persist)) |>
    group_by(horizon) |>
    summarise(cor_persistence = cor(mu_med, persist),
              cor_truth = cor(mu_med, .truth), .groups = "drop")
  cat("\n--- is this persistence in disguise? ---\n")
  print(as.data.frame(lag_chk), row.names = FALSE, digits = 3)
  if (any(lag_chk$cor_persistence > 0.95))
    cat("  correlation with persistence above 0.95: the model is reproducing\n",
        "  the last observed value. A random walk cannot do otherwise past the\n",
        "  origin; try interaction = 'harmonic_ar1' with add_trend = FALSE.\n",
        sep = "")

  cat("\nskill > 0 beats the baseline. Persistence is hard to beat at one week;\n",
      "SEASONAL NAIVE is the meaningful test for a seasonal disease.\n",
      "Check coverage before MAE: intervals that are not calibrated make the\n",
      "point accuracy hard to interpret.\n", sep = "")

  by_state <- d |> filter(!cen_t, is.finite(.truth)) |>
    group_by(state, horizon) |>
    summarise(n = n(), MAE = mean(abs(mu_med - .truth)),
              coverage = mean(.truth >= mu_lo & .truth <= mu_hi),
              .groups = "drop")
  invisible(list(overall = out, by_state = by_state, data = d))
}


## ---------------------------------------------------------------------------
## figures
## ---------------------------------------------------------------------------
plot_nowcast_inla <- function(nc, dat, horizon = 1, states = NULL,
                              file = NULL, base = 14, ncol = 3) {
  d <- nc |> filter(horizon == !!horizon)
  if (!is.null(states)) d <- filter(d, state %in% states)
  obs <- dat$data |>
    filter(state %in% unique(d$state), !censored, is.finite(y_raw),
           week >= min(d$week) - 180, week <= max(d$week))
  fl <- function(x) pmax(x, 0.5)

  p <- ggplot() +
    geom_point(data = obs, aes(week, fl(y_raw)), colour = "black", size = .7) +
    geom_ribbon(data = d, aes(week, ymin = fl(mu_lo), ymax = fl(mu_hi),
                              group = state), fill = "#2166ac", alpha = .25) +
    geom_line(data = d, aes(week, fl(mu_med), group = state),
              colour = "#2166ac", linewidth = .8) +
    scale_y_log10() +
    facet_wrap(~ state, ncol = ncol, scales = "free_y") +
    labs(title = paste0(horizon, "-week-ahead nowcast"),
         subtitle = paste("blue = predictive median and 95% interval, one point",
                          "per origin;\nblack = observed counts"),
         x = NULL, y = "weekly RSV encounters (log scale)") +
    theme_minimal(base_size = base) +
    theme(strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base + 3))
  if (!is.null(file)) { ggsave(file, p, width = 13, height = 8, dpi = 150)
                        message("wrote ", file) }
  p
}


#' Skill and coverage against horizon, both baselines shown.
plot_nowcast_scores <- function(sc, file = NULL, base = 14) {
  o <- sc$overall
  a <- o |> select(horizon, MAE, MAE_persistence, MAE_seasonal) |>
    pivot_longer(-horizon, names_to = "model", values_to = "MAE") |>
    mutate(model = recode(model, MAE = "model",
                          MAE_persistence = "persistence",
                          MAE_seasonal = "seasonal naive"))
  p1 <- ggplot(a, aes(factor(horizon), MAE, fill = model)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c(model = "#2166ac", persistence = "grey55",
                                 `seasonal naive` = "#b2182b"), name = NULL) +
    labs(title = "Nowcast error by horizon",
         subtitle = "the model must beat BOTH baselines to claim skill",
         x = "weeks ahead", y = "MAE (counts)") +
    theme_minimal(base_size = base) + theme(legend.position = "top")

  p2 <- ggplot(o, aes(factor(horizon), coverage)) +
    geom_col(fill = "#4c72b0") +
    geom_hline(yintercept = .95, linetype = "dashed", colour = "red") +
    ylim(0, 1) +
    labs(title = "95% predictive interval coverage",
         subtitle = "dashed = nominal", x = "weeks ahead", y = "coverage") +
    theme_minimal(base_size = base)

  if (!is.null(file)) {
    ggsave(sub("\\.png$", "_mae.png", file), p1, width = 7, height = 4.5, dpi = 150)
    ggsave(sub("\\.png$", "_coverage.png", file), p2, width = 7, height = 4.5, dpi = 150)
    message("wrote two score figures")
  }
  list(mae = p1, coverage = p2)
}


## ---------------------------------------------------------------------------
## every horizon on one panel
## ---------------------------------------------------------------------------
#' All horizons together, shaded by how far ahead the forecast is.
#'
#' Four separate figures make the reader hold one in their head while looking at
#' the next. On one panel the degradation with horizon is immediate: the lines
#' spread apart and the bands widen as h grows, and where they diverge tells you
#' at which point in the season the extra weeks start to cost something.
#'
#' Ribbons for four horizons overlaid become unreadable, so by default only the
#' shortest and longest are shaded -- they bracket the rest. `ribbons = "all"`
#' draws every band at low opacity, `ribbons = "none"` drops them entirely and
#' is the clearest option when the point is the medians.
plot_nowcast_horizons <- function(nc, dat, states = NULL, file = NULL,
                                  base = 14, ncol = 3,
                                  ribbons = c("extremes", "all", "none"),
                                  width = 13, height = 8) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
  ribbons <- match.arg(ribbons)

  d <- nc
  if (!is.null(states)) d <- filter(d, state %in% states)
  d <- d |> mutate(h = factor(horizon, levels = sort(unique(horizon))))
  hs <- levels(d$h)

  obs <- dat$data |>
    filter(state %in% unique(d$state), !censored, is.finite(y_raw),
           week >= min(d$week) - 220, week <= max(d$week))
  fl <- function(x) pmax(x, 0.5)

  ## sequential blues: light for one week ahead, dark for the longest
  pal <- setNames(colorRampPalette(c("#9ecae1", "#08306b"))(length(hs)), hs)

  rib <- switch(ribbons,
                extremes = hs[c(1, length(hs))],
                all      = hs,
                none     = character(0))

  p <- ggplot()
  ## widest band first, so the narrow ones stay visible on top
  for (hh in rev(rib)) {
    dh <- filter(d, h == hh)
    p <- p + geom_ribbon(data = dh,
                         aes(week, ymin = fl(mu_lo), ymax = fl(mu_hi),
                             group = interaction(state, h), fill = h),
                         alpha = if (ribbons == "all") .12 else .18)
  }
  p <- p +
    geom_point(data = obs, aes(week, fl(y_raw)), colour = "grey15", size = .7) +
    geom_line(data = d, aes(week, fl(mu_med), colour = h,
                            group = interaction(state, h)), linewidth = .6) +
    scale_colour_manual(values = pal, name = "weeks ahead") +
    scale_fill_manual(values = pal, name = "weeks ahead", guide = "none") +
    scale_y_log10() +
    facet_wrap(~ state, ncol = ncol, scales = "free_y") +
    labs(title = "Nowcast at every horizon",
         subtitle = paste0("black = observed;  lines = predictive median, ",
                           "light to dark with increasing horizon",
                           if (ribbons == "extremes")
                             ";\nshading = 95% interval at the shortest and "
                           else if (ribbons == "all") ";\nshading = 95% intervals"
                           else "",
                           if (ribbons == "extremes") "longest horizon" else ""),
         x = NULL, y = "weekly RSV encounters (log scale)") +
    theme_minimal(base_size = base) +
    theme(strip.text = element_text(face = "bold"), legend.position = "top")

  if (!is.null(file)) { ggsave(file, p, width = width, height = height, dpi = 300)
                        message("wrote ", file) }
  p
}


#' The same idea for one state, large enough to read the detail.
plot_nowcast_one_state <- function(nc, dat, state_name, file = NULL, base = 14) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
  d <- nc |> filter(state == state_name) |>
    mutate(h = factor(horizon, levels = sort(unique(horizon))))
  obs <- dat$data |> filter(state == state_name, !censored, is.finite(y_raw),
                            week >= min(d$week) - 220, week <= max(d$week))
  hs <- levels(d$h)
  pal <- setNames(colorRampPalette(c("#9ecae1", "#08306b"))(length(hs)), hs)
  fl <- function(x) pmax(x, 0.5)

  p <- ggplot() +
    geom_ribbon(data = filter(d, h == tail(hs, 1)),
                aes(week, ymin = fl(mu_lo), ymax = fl(mu_hi)),
                fill = tail(pal, 1), alpha = .16) +
    geom_ribbon(data = filter(d, h == hs[1]),
                aes(week, ymin = fl(mu_lo), ymax = fl(mu_hi)),
                fill = pal[1], alpha = .22) +
    geom_point(data = obs, aes(week, fl(y_raw)), colour = "grey15", size = 1.1) +
    geom_line(data = d, aes(week, fl(mu_med), colour = h), linewidth = .8) +
    scale_colour_manual(values = pal, name = "weeks ahead") +
    scale_y_log10() +
    labs(title = paste0("Nowcast at every horizon: ", state_name),
         subtitle = paste("black = observed;  lines light to dark with horizon;",
                          "shading = 95% interval at 1 and at the longest"),
         x = NULL, y = "weekly RSV encounters (log scale)") +
    theme_minimal(base_size = base) + theme(legend.position = "top")
  if (!is.null(file)) { ggsave(file, p, width = 9, height = 5.5, dpi = 300)
                        message("wrote ", file) }
  p
}


## ---------------------------------------------------------------------------
## per-state scores
## ---------------------------------------------------------------------------
#' Score one structure state by state, at each horizon.
#'
#' The pooled MAE is dominated by the largest states: Texas contributes errors
#' in the hundreds while Vermont contributes errors in single digits, so a
#' national figure mostly describes how well the model does in Texas. Per-state
#' scores show where the forecast is usable, and relative error makes states of
#' very different size comparable.
#'
#' Each state's own baselines are computed too, because "beats persistence" is a
#' claim that can hold nationally while failing in half the states -- and it is
#' the per-state version a health department would care about.
#'
#' SCORING THE DISTRIBUTION, NOT JUST THE POINT
#'
#' MAE judges the median and ignores the interval, so a model can win on MAE
#' while being badly calibrated. Two proper scores are added:
#'
#'   WIS  the weighted interval score. CRPS cannot be computed here because the
#'        full predictive distribution is not stored -- only the median and a
#'        95% interval -- and WIS is the quantile-based approximation to CRPS
#'        used by the COVID and influenza Forecast Hubs. With one central
#'        interval at level 1-alpha,
#'
#'          IS   = (u - l) + (2/alpha)(l - y)1{y<l} + (2/alpha)(y - u)1{y>u}
#'          WIS  = [ |y - m|/2 + (alpha/2) IS ] / 1.5
#'
#'        It rewards a narrow interval and penalises one that misses, so a model
#'        cannot buy coverage by being vague. Lower is better, and it is on the
#'        scale of the counts, like MAE.
#'
#'   CRPS_pois  exact CRPS under Poisson(mu_med). This IS a closed-form CRPS,
#'        but for a different object: it treats the fitted mean as if it were
#'        known, so it ignores the parameter uncertainty that the interval
#'        represents. Reported as a cross-check, not as the headline.
#'
#' @param nc      one structure's predictions, from rolling_nowcast()
#' @param dat     stcar_prepare() output, for the baselines
#' @param cen_hi  the censoring ceiling, for the suppressed-origin fallback
#' @param alpha   the interval's nominal miss rate (0.05 for a 95% interval)
score_nowcast_by_state <- function(nc, dat, cen_hi = 10, min_n = 5,
                                   alpha = 0.05,
                                   out_dir = NULL, base = 13) {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(tidyr)
  })

  src <- dat$data |> transmute(state, week, y = as.numeric(y_raw), censored)
  pers <- src |> transmute(state, origin = week,
                           persist = ifelse(censored, cen_hi / 2, y))
  sn <- src |> transmute(state, w52 = week,
                         seas = ifelse(censored, cen_hi / 2, y))

  d <- nc |>
    left_join(src |> transmute(state, week, cen_t = censored),
              by = c("state", "week")) |>
    left_join(pers, by = c("state", "origin")) |>
    mutate(w52 = week - 7 * 52) |>
    left_join(sn, by = c("state", "w52")) |>
    filter(!cen_t, is.finite(.truth))

  d <- d |>
    mutate(WIS  = .wis1(.truth, mu_med, mu_lo, mu_hi, alpha),
           CRPS_pois = .crps_pois(.truth, mu_med),
           ## the same score for the baselines, which as point forecasts have
           ## no interval: WIS then reduces to |y - m|, so the comparison is
           ## fair only in that it scores them on what they actually provide
           WIS_persistence = abs(.truth - persist),
           WIS_seasonal    = abs(.truth - seas))

  by_state <- d |>
    group_by(state, horizon) |>
    summarise(
      n        = n(),
      mean_obs = mean(.truth),
      MAE      = mean(abs(mu_med - .truth)),
      RMSE     = sqrt(mean((mu_med - .truth)^2)),
      WIS      = mean(WIS),
      CRPS_pois = mean(CRPS_pois),
      ## relative error is what makes a small state comparable with a large one
      rel_MAE  = mean(abs(mu_med - .truth) / pmax(.truth, 1)),
      bias     = mean(mu_med - .truth),
      corr     = if (n() > 2 && sd(.truth) > 0 && sd(mu_med) > 0)
                   cor(mu_med, .truth) else NA_real_,
      coverage = mean(.truth >= mu_lo & .truth <= mu_hi),
      MAE_persistence = mean(abs(persist - .truth), na.rm = TRUE),
      MAE_seasonal    = mean(abs(seas - .truth), na.rm = TRUE),
      .groups = "drop") |>
    mutate(skill_vs_persistence = 1 - MAE / MAE_persistence,
           skill_vs_seasonal    = 1 - MAE / MAE_seasonal,
           ## skill on the proper score, which is the one to lead with: it
           ## credits a sharp well-placed interval and penalises a vague one,
           ## where MAE is blind to both
           wis_skill_persistence = 1 - WIS / MAE_persistence)

  thin <- by_state |> filter(n < min_n) |> distinct(state) |> pull(state)
  if (length(thin))
    message("fewer than ", min_n, " scorable targets in: ",
            paste(thin, collapse = ", "),
            "\n  their scores are unstable and are flagged in the output")
  by_state$sparse <- by_state$n < min_n

  ## ---- what the pooled number hides -------------------------------------
  cat("\n=========== PER-STATE SCORES ===========\n")
  summ <- by_state |> filter(!sparse) |>
    group_by(horizon) |>
    summarise(states = n(),
              median_MAE = median(MAE),
              median_WIS = median(WIS),
              median_CRPS_pois = median(CRPS_pois),
              median_rel_MAE = median(rel_MAE),
              median_coverage = median(coverage),
              beat_persistence = mean(skill_vs_persistence > 0),
              beat_seasonal = mean(skill_vs_seasonal > 0),
              beat_persistence_WIS = mean(wis_skill_persistence > 0),
              .groups = "drop")
  print(as.data.frame(summ), row.names = FALSE, digits = 3)
  cat("\nbeat_persistence is the SHARE OF STATES where the model wins, not the\n",
      "pooled comparison. A model can beat persistence nationally while losing\n",
      "in half the states, and it is the per-state version that matters to a\n",
      "health department.\n", sep = "")

  for (h in sort(unique(by_state$horizon))) {
    b <- by_state |> filter(horizon == h, !sparse) |> arrange(desc(rel_MAE))
    cat("\nhorizon ", h, " -- worst five by relative error:\n", sep = "")
    print(b |> select(state, n, mean_obs, MAE, WIS, rel_MAE, coverage,
                      skill_vs_persistence) |> head(5) |> as.data.frame(),
          row.names = FALSE, digits = 3)
  }

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(by_state, file.path(out_dir, "nowcast_scores_by_state.csv"))
    readr::write_csv(summ, file.path(out_dir, "nowcast_scores_by_state_summary.csv"))

    ## relative error by state, one bar per horizon
    p1 <- by_state |> filter(!sparse) |>
      ggplot(aes(reorder(state, rel_MAE), rel_MAE, fill = factor(horizon))) +
      geom_col(position = "dodge") + coord_flip() +
      scale_fill_manual(values = colorRampPalette(c("#9ecae1", "#08306b"))(
        length(unique(by_state$horizon))), name = "weeks ahead") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(title = "Relative nowcast error by state",
           subtitle = "mean |predicted - observed| / observed",
           x = NULL, y = NULL) +
      theme_minimal(base_size = base) +
      theme(axis.text.y = element_text(size = 6), legend.position = "top")
    ggsave(file.path(out_dir, "nowcast_rel_error_by_state.png"), p1,
           width = 9, height = 11, dpi = 300)

    ## where does the model beat persistence?
    p2 <- by_state |> filter(!sparse) |>
      ggplot(aes(reorder(state, skill_vs_persistence), skill_vs_persistence,
                 fill = factor(horizon))) +
      geom_col(position = "dodge") + coord_flip() +
      geom_hline(yintercept = 0, linewidth = .4) +
      scale_fill_manual(values = colorRampPalette(c("#9ecae1", "#08306b"))(
        length(unique(by_state$horizon))), name = "weeks ahead") +
      labs(title = "Skill against persistence, by state",
           subtitle = "above zero the model beats carrying the last count forward",
           x = NULL, y = "skill") +
      theme_minimal(base_size = base) +
      theme(axis.text.y = element_text(size = 6), legend.position = "top")
    ggsave(file.path(out_dir, "nowcast_skill_by_state.png"), p2,
           width = 9, height = 11, dpi = 300)

    ## WIS by state: the proper score, so this is the panel to lead with
    p4 <- by_state |> filter(!sparse) |>
      ggplot(aes(reorder(state, WIS), WIS, fill = factor(horizon))) +
      geom_col(position = "dodge") + coord_flip() +
      scale_fill_manual(values = colorRampPalette(c("#9ecae1", "#08306b"))(
        length(unique(by_state$horizon))), name = "weeks ahead") +
      labs(title = "Weighted interval score by state",
           subtitle = paste("a proper score: rewards a narrow interval,",
                            "penalises one that misses.  Lower is better."),
           x = NULL, y = "WIS (counts)") +
      theme_minimal(base_size = base) +
      theme(axis.text.y = element_text(size = 6), legend.position = "top")
    ggsave(file.path(out_dir, "nowcast_wis_by_state.png"), p4,
           width = 9, height = 11, dpi = 300)

    ## coverage, which should sit near 0.95 everywhere
    p3 <- by_state |> filter(!sparse) |>
      ggplot(aes(factor(horizon), coverage)) +
      geom_hline(yintercept = .95, linetype = "dashed", colour = "#b2182b") +
      geom_boxplot(fill = "#4c78a8", alpha = .3, outlier.size = .7) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "95% interval coverage across states",
           subtitle = "dashed = nominal;  each point is one state",
           x = "weeks ahead", y = "coverage") +
      theme_minimal(base_size = base)
    ggsave(file.path(out_dir, "nowcast_coverage_by_state.png"), p3,
           width = 7, height = 4.5, dpi = 300)

    message("wrote per-state scores and three figures to ", out_dir)
  }
  invisible(list(by_state = by_state, summary = summ, data = d))
}


#' Per-state scores for several structures at once, matched on origins.
score_by_state_all <- function(nc_all, dat, cen_hi = 10, out_dir = NULL) {
  bind_rows(lapply(names(nc_all), function(ix) {
    s <- score_nowcast_by_state(nc_all[[ix]], dat, cen_hi = cen_hi,
                                out_dir = NULL)
    s$by_state |> mutate(interaction = ix)
  })) -> tab
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(tab, file.path(out_dir, "nowcast_by_state_all_structures.csv"))
  }
  cat("\n=========== BEST STRUCTURE PER STATE (horizon 4) ===========\n")
  best <- tab |> filter(horizon == 4, !sparse) |>
    group_by(state) |> slice_min(MAE, n = 1, with_ties = FALSE) |>
    ungroup() |> count(interaction, sort = TRUE)
  print(as.data.frame(best), row.names = FALSE)
  cat("\nIf one structure wins nearly everywhere the national choice is safe;\n",
      "if the wins are split, the best model depends on the state and a single\n",
      "national ranking hides that.\n", sep = "")
  invisible(tab)
}


## ---------------------------------------------------------------------------
## comparing structures on a chosen score, in chosen states
## ---------------------------------------------------------------------------
#' Bar plot of any per-state score, across structures and horizons.
#'
#' `score_by_state_all()` produces a table with one row per structure, state and
#' horizon. Which number matters depends on the question, so the metric is an
#' argument rather than fixed:
#'
#'   WIS       the proper score. Credits a sharp interval that contains the
#'             truth and penalises both vagueness and overconfidence, so this is
#'             the one to lead with when the intervals are part of the claim.
#'   MAE       the median only. Familiar, on the scale of counts, but blind to
#'             whether the interval was any good.
#'   rel_MAE   MAE divided by the observed count, which is what makes a small
#'             state comparable with a large one. Use it whenever several states
#'             of different size are on the same axis -- otherwise the panel is
#'             really a plot of state size.
#'   coverage  should sit near 0.95; far from it, the other scores are hard to
#'             read.
#'   skill_vs_persistence, wis_skill_persistence
#'             above zero the model beats carrying the last count forward.
#'
#' Lower is better for WIS, MAE, RMSE and rel_MAE; higher for the skills and
#' for coverage up to its nominal level. The function orients the reference line
#' and the caption accordingly.
#'
#' @param tab     score_by_state_all() output
#' @param metric  which column to plot
#' @param states  NULL for every state pooled, or a character vector to facet
#' @param horizons  which horizons to show
#' @param structures  which structures, and in what order; NULL for all
plot_score_comparison <- function(tab, metric = "WIS", states = NULL,
                                  horizons = NULL, structures = NULL,
                                  file = NULL, base = 13, ncol = 3,
                                  drop_sparse = TRUE) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

  if (!metric %in% names(tab))
    stop("`metric` must be one of: ",
         paste(setdiff(names(tab), c("state", "horizon", "interaction",
                                     "sparse", "n")), collapse = ", "),
         call. = FALSE)

  d <- tab
  if (drop_sparse && "sparse" %in% names(d)) d <- filter(d, !sparse)
  if (!is.null(states))     d <- filter(d, state %in% states)
  if (!is.null(horizons))   d <- filter(d, horizon %in% horizons)
  if (!is.null(structures)) d <- filter(d, interaction %in% structures)
  if (!nrow(d)) stop("nothing left after filtering", call. = FALSE)

  ix <- if (is.null(structures)) unique(d$interaction) else structures
  d$interaction <- factor(d$interaction, levels = ix)
  pal <- setNames(colorRampPalette(c("#9ecae1", "#08306b"))(length(ix)), ix)

  lower_better <- metric %in% c("WIS", "MAE", "RMSE", "rel_MAE", "CRPS_pois",
                                "MAE_persistence", "MAE_seasonal")
  is_skill <- grepl("skill", metric)

  ## With no states given, the states are POOLED rather than averaged: a mean of
  ## per-state MAEs would weight Vermont equally with Texas, which is a
  ## different quantity from the national error and is easy to misread as it.
  if (is.null(states)) {
    d <- d |> group_by(interaction, horizon) |>
      summarise(value = if (metric %in% c("MAE", "WIS", "RMSE", "CRPS_pois"))
                  stats::weighted.mean(.data[[metric]], n)
                else median(.data[[metric]]),
                .groups = "drop")
    sub <- if (metric %in% c("MAE", "WIS", "RMSE", "CRPS_pois"))
      "pooled across states, weighted by the number of scored targets"
    else "median across states"
    fac <- NULL
  } else {
    d <- d |> transmute(interaction, horizon, state,
                        value = .data[[metric]])
    sub <- paste("one panel per state;", length(unique(d$state)), "states")
    fac <- facet_wrap(~ state, ncol = ncol, scales = "free_y")
  }

  p <- ggplot(d, aes(factor(horizon), value, fill = interaction)) +
    geom_col(position = position_dodge(width = .8), width = .75) +
    scale_fill_manual(values = pal, name = NULL) +
    labs(title = paste0(metric, " by temporal structure"),
         subtitle = paste0(sub, ";  ",
                           if (is_skill) "above zero beats the baseline"
                           else if (lower_better) "lower is better"
                           else "higher is better"),
         x = "weeks ahead", y = metric) +
    theme_minimal(base_size = base) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold"))

  if (is_skill) p <- p + geom_hline(yintercept = 0, linewidth = .4)
  if (metric == "coverage")
    p <- p + geom_hline(yintercept = .95, linetype = "dashed",
                        colour = "#b2182b") +
      coord_cartesian(ylim = c(0, 1))
  if (!is.null(fac)) p <- p + fac

  if (!is.null(file)) {
    h <- if (is.null(states)) 5 else max(5, 2.6 * ceiling(length(unique(d$state)) / ncol))
    ggsave(file, p, width = if (is.null(states)) 8.5 else 4.2 * ncol,
           height = h, dpi = 300)
    message("wrote ", file)
  }

  ## which structure wins, so the ranking does not have to be read off bars
  best <- d |> group_by(across(any_of(c("state", "horizon")))) |>
    slice(if (lower_better) which.min(value) else which.max(value)) |>
    ungroup()
  cat("\nbest structure on ", metric, ":\n", sep = "")
  print(best |> count(interaction, sort = TRUE) |> as.data.frame(),
        row.names = FALSE)

  invisible(list(plot = p, data = d, best = best))
}


## ---------------------------------------------------------------------------
## one fan per origin
## ---------------------------------------------------------------------------
#' Each origin's own four-week trajectory, drawn as a short fan.
#'
#' `plot_nowcast_horizons()` joins all the h=1 predictions into one line, all
#' the h=2 predictions into another, and so on. That shows how accuracy decays
#' with horizon, but it hides what a user actually sees: at any given week, the
#' forecast is a short trajectory running four weeks forward, and the next week
#' brings a new one.
#'
#' This draws that instead -- one line and one band per origin, spanning
#' horizons 1 to 4. Successive fans overlap, and where they disagree the model
#' changed its mind as data arrived. A fan that repeatedly starts above the
#' observations and bends down is a model being pulled back each week; one that
#' tracks smoothly is a model whose view is stable.
#'
#' @param origin_every  draw every k-th origin. Overlapping fans are the point,
#'   but with weekly origins and four-week horizons each week carries four, and
#'   the panel becomes solid. Thinning is a readability choice and should be
#'   stated.
#' @param probs  extra interval levels, if the predictions carry them. With only
#'   a 95% interval stored, one band is all that can be drawn; the nested bands
#'   in forecast-hub figures need the intermediate quantiles to be saved at
#'   prediction time.
#' @param origins  an explicit set of origin dates to draw, or a length-2
#'   vector giving a range. Overrides `origin_every`. Useful for showing one
#'   part of the season at full density -- the epidemic rise, say -- rather
#'   than thinning the whole year uniformly.
plot_nowcast_fans <- function(nc, dat, states = NULL, file = NULL,
                              base = 12, ncol = 6, origin_every = 1,
                              origins = NULL,
                              colour = "#2b4a9b", linear = TRUE,
                              width = 16, height = 10) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

  d <- nc
  if (!is.null(states)) d <- filter(d, state %in% states)
  if (!nrow(d)) stop("nothing to plot", call. = FALSE)

  ## Thinning is a READABILITY choice and changes nothing that is scored: the
  ## scores were computed before this function was called and use every origin.
  ## It only decides how many fans are drawn.
  all_og <- sort(unique(d$origin))
  if (!is.null(origins)) {
    keep <- if (length(origins) == 2 && inherits(as.Date(origins), "Date"))
      all_og[all_og >= as.Date(origins)[1] & all_og <= as.Date(origins)[2]]
    else all_og[all_og %in% as.Date(origins)]
    if (!length(keep)) stop("no origins in the range given", call. = FALSE)
    d <- filter(d, origin %in% keep)
    message("drawing ", length(keep), " of ", length(all_og), " origins (",
            format(min(keep)), " to ", format(max(keep)), ")")
  } else if (origin_every > 1) {
    keep <- all_og[seq(1, length(all_og), by = origin_every)]
    d <- filter(d, origin %in% keep)
    message("drawing every ", origin_every, "th origin (", length(keep),
            " of ", length(all_og), ")")
  }

  obs <- dat$data |>
    filter(state %in% unique(d$state), !censored, is.finite(y_raw),
           week >= min(d$week), week <= max(d$week))

  ## A linear axis matches the forecast-hub convention and makes the fans
  ## legible; a log axis would compress the epidemic peak, which is where the
  ## disagreement between fans is most informative. Small states then have
  ## little vertical range, which is why the scales are free.
  p <- ggplot() +
    geom_ribbon(data = d, aes(week, ymin = mu_lo, ymax = mu_hi, group = origin),
                fill = colour, alpha = .13) +
    geom_line(data = d, aes(week, mu_med, group = origin),
              colour = colour, linewidth = .4) +
    geom_point(data = d, aes(week, mu_med, group = origin),
               colour = colour, size = .35) +
    geom_line(data = obs, aes(week, y_raw), colour = "black", linewidth = .45) +
    geom_point(data = obs, aes(week, y_raw), colour = "black", size = .5) +
    facet_wrap(~ state, ncol = ncol, scales = "free_y") +
    labs(title = "Forecasts as they were made",
         subtitle = paste("black = observed;  each blue fan is one origin's",
                          "trajectory over the next four weeks,\nwith its 95%",
                          "interval. Overlapping fans show the model revising",
                          "as data arrives."),
         x = NULL, y = "weekly RSV encounters") +
    theme_minimal(base_size = base) +
    theme(strip.text = element_text(face = "bold", size = base - 2),
          axis.text = element_text(size = base - 4),
          panel.grid.minor = element_blank())
  if (!linear) p <- p + scale_y_log10()

  if (!is.null(file)) { ggsave(file, p, width = width, height = height, dpi = 300)
                        message("wrote ", file) }
  p
}


#' The same for one state, large enough to see individual fans.
plot_nowcast_fans_one <- function(nc, dat, state_name, file = NULL,
                                  base = 14, origin_every = 1, origins = NULL,
                                  colour = "#2b4a9b", linear = TRUE) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
  d <- nc |> filter(state == state_name)
  all_og <- sort(unique(d$origin))
  if (!is.null(origins)) {
    keep <- if (length(origins) == 2)
      all_og[all_og >= as.Date(origins)[1] & all_og <= as.Date(origins)[2]]
    else all_og[all_og %in% as.Date(origins)]
    d <- filter(d, origin %in% keep)
  } else if (origin_every > 1) {
    d <- filter(d, origin %in% all_og[seq(1, length(all_og), by = origin_every)])
  }
  obs <- dat$data |> filter(state == state_name, !censored, is.finite(y_raw),
                            week >= min(d$week), week <= max(d$week))
  p <- ggplot() +
    geom_ribbon(data = d, aes(week, ymin = mu_lo, ymax = mu_hi, group = origin),
                fill = colour, alpha = .15) +
    geom_line(data = d, aes(week, mu_med, group = origin),
              colour = colour, linewidth = .6) +
    geom_point(data = d, aes(week, mu_med), colour = colour, size = .9) +
    geom_line(data = obs, aes(week, y_raw), colour = "black", linewidth = .7) +
    geom_point(data = obs, aes(week, y_raw), colour = "black", size = 1.4) +
    labs(title = paste0("Forecasts as they were made: ", state_name),
         subtitle = paste("each fan is one origin's next four weeks, with its",
                          "95% interval"),
         x = NULL, y = "weekly RSV encounters") +
    theme_minimal(base_size = base)
  if (!linear) p <- p + scale_y_log10()
  if (!is.null(file)) { ggsave(file, p, width = 10, height = 5.5, dpi = 300)
                        message("wrote ", file) }
  p
}
