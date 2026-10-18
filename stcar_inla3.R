
suppressPackageStartupMessages({ library(dplyr) })

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Which censored likelihoods does this INLA build provide?
inla_censored_families <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) stop("INLA not installed")
  fams <- try(INLA::inla.list.models("likelihood"), silent = TRUE)
  if (inherits(fams, "try-error")) return(character(0))
  nm <- if (is.list(fams)) names(fams) else as.character(fams)
  grep("cen", nm, value = TRUE, ignore.case = TRUE)
}


#' @param dat          stcar_prepare() output
#' @param climate      time-varying covariates z_it
#' @param use_pop      include the time-invariant x_i (see warning 2)
#' @param slope_model  "iid" (independent state slopes) or "besag" (spatially
#'                     correlated slopes); "none" for a common slope only
#' @param family       "nbinomial" (censored if this build supports it),
#'                     "cenpoisson", or "poisson"
#' @param interaction  "type2", "type4", or "none"
stcar_fit_inla3 <- function(dat,
                            climate = c("abs_humidity_lag1", "temp_lag1"),
                            use_pop = TRUE,
                            mobility = c("none", "anomaly", "raw"),
                            mob_lag = 1,
                            extra_fixed = NULL,
                            hyper_priors = NULL,
                            quantiles = NULL,
                            harmonics = FALSE,
                            slope_model = c("iid", "besag", "none"),
                            ## Zero-inflated families matter when a sizeable
                            ## share of the panel is EXACTLY zero. A log-link
                            ## Poisson has mu > 0 everywhere, so its lower
                            ## predictive quantile is strictly positive and an
                            ## observed zero is missed by construction, however
                            ## good the mean is. In the CDC hub's NHSN panel
                            ## 26% of state-weeks are zero -- concentrated in
                            ## Vermont, Wyoming, DC, Montana and the Dakotas --
                            ## and those alone account for the whole gap
                            ## between 83% and 95% coverage.
                            ##
                            ## Type 1 ("...poisson1") models the zero
                            ## probability as a free parameter: some weeks have
                            ## no RSV admissions at all, rather than a small
                            ## positive rate that happened to round down.
                            family = c("nbinomial", "cenpoisson", "poisson",
                                       "zeroinflatedpoisson1",
                                       "zeroinflatednbinomial1"),
                            interaction = c("type2", "type4", "harmonic",
                                            "harmonic_ar1", "harmonic_season",
                                            "monthly", "ar1", "none"),
                            censored = TRUE, subset_states = NULL,
                            add_trend = TRUE, light = FALSE,
                            cen_lo = 1, cen_hi = 10, verbose = FALSE) {

  slope_model <- match.arg(slope_model)
  mobility    <- match.arg(mobility)
  family      <- match.arg(family)
  interaction <- match.arg(interaction)
  if (!requireNamespace("INLA", quietly = TRUE)) stop("INLA not installed")

  d <- dat$data; W <- dat$W
  if (!is.null(subset_states)) {
    d <- dplyr::filter(d, state %in% subset_states)
    W <- W[subset_states, subset_states, drop = FALSE]
    states <- subset_states
  } else states <- dat$states

  d <- d |>
    mutate(i = match(state, states),
           t = match(week, sort(unique(week))),
           i2 = i, t2 = t)

  ## ---- response and family --------------------------------------------------
  cenfams <- inla_censored_families()
  fam <- family; cf <- list()
  if (censored) {
    if (family == "nbinomial") {
      if (any(grepl("nbinomial", cenfams))) {
        fam <- grep("nbinomial", cenfams, value = TRUE)[1]
        cf  <- list(cenpoisson.I = c(cen_lo, cen_hi))   # name may differ; see below
        message("using censored NB family: ", fam)
        message("  VERIFY the interval argument name with inla.doc('", fam, "')")
      } else {
        warning("this INLA build has no censored negative binomial (found: ",
                if (length(cenfams)) paste(cenfams, collapse = ", ") else "none",
                ").\n  Falling back to cenpoisson. Note that with phi_it in the",
                " model the interaction already carries the overdispersion, so",
                " little is lost.", call. = FALSE)
        fam <- "cenpoisson"; cf <- list(cenpoisson.I = c(cen_lo, cen_hi))
      }
    } else if (family == "cenpoisson") {
      cf <- list(cenpoisson.I = c(cen_lo, cen_hi))
    } else if (grepl("^zeroinflated", family)) {
      ## Censoring and zero inflation answer different questions -- "the count
      ## is somewhere in an interval" against "the count is exactly zero" --
      ## and INLA has no family that does both. Censoring wins here because it
      ## is the harder constraint to approximate.
      warning("no censored zero-inflated family exists; censoring takes ",
              "precedence and the zero inflation is dropped.\n  Use ",
              "censored = FALSE with a zero-inflated family on an ",
              "uncensored panel.", call. = FALSE)
      fam <- "cenpoisson"; cf <- list(cenpoisson.I = c(cen_lo, cen_hi))
    }
    d$y_inla <- as.integer(ifelse(d$censored, cen_lo, d$y_raw))
  } else {
    d$y_inla <- as.integer(ifelse(d$censored, NA, d$y_raw))
    fam <- sub("^cen", "", fam)
  }

  ## ---- covariates -----------------------------------------------------------
  zvars <- climate[climate %in% names(d)]
  if (length(zvars) < length(climate)) {
    miss <- setdiff(climate, zvars)
    warning("climate covariate(s) NOT FOUND and therefore EXCLUDED: ",
            paste(miss, collapse = ", "),
            "\n  Columns are named <variable>_lag<k> and are built by ",
            "stcar_prepare(climate_lags = ...).",
            "\n  Available here: ",
            paste(grep("_lag[0-9]+$", names(d), value = TRUE), collapse = ", "),
            call. = FALSE)
  }

  scale1 <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x[!is.finite(x)] <- NA
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    s <- stats::sd(x); if (s > 0) (x - mean(x)) / s else rep(0, length(x))
  }
  for (v in zvars) d[[v]] <- scale1(d[[v]])

  fixed <- character(0)
  if (use_pop) {
    if ("log_pop" %in% names(d)) {
      d$x_pop <- scale1(d$log_pop); fixed <- c(fixed, "x_pop")
      message("including x_i = log(population). NOTE: constant within state and",
              " therefore confounded with b_i; treat beta_pop with caution.")
    } else message("log_pop not available; skipping x_i")
  }
  ## ---- MOBILITY -------------------------------------------------------------
  ## "raw"     mob_pressure_lag1: sum_j M_ij * (RSV rate in state j last week).
  ##           CAUTION: built FROM other states' RSV, so it shares the national
  ##           seasonality of the target state. Its raw correlation with the
  ##           local rate was 0.76, almost all of it common seasonal timing, and
  ##           in the NIMBLE fit it did not converge (R-hat 3.7, ESS 8).
  ##
  ## "anomaly" mob_anom_lag1 = mob_pressure_lag1 minus the national mean that
  ##           week, i.e. imported pressure RELATIVE to the country as a whole.
  ##           The seasonal component is removed, so what remains is closer to a
  ##           genuine spatial spillover. Correlation fell to 0.32.
  ##           Computed here if the column is absent.
  ##
  ## Either way the term competes with gamma_t and phi_it, which already carry
  ## shared temporal signal; treat a large coefficient sceptically.
  ## `mob_lag` is how many weeks back the other-state signal is taken from.
  ## For FITTING, 1 is natural. For FORECASTING it must be at least the horizon,
  ## or the covariate on a future row is built from counts that have not
  ## happened yet at the origin -- which leaks and would flatter the forecast.
  if (mobility != "none") {
    pv <- paste0("mob_pressure_lag", mob_lag)
    av <- paste0("mob_anom_lag", mob_lag)
    mv <- if (mobility == "anomaly") av else pv
    if (!pv %in% names(d))
      stop(pv, " is not in the data. Rebuild with ",
           "stcar_prepare(mob_lags = 1:", max(mob_lag, 8), ").", call. = FALSE)
    if (mobility == "anomaly" && !(av %in% names(d))) {
      d[[av]] <- d[[pv]]
      d <- d |> group_by(week) |>
        mutate(!!av := .data[[av]] - mean(.data[[av]], na.rm = TRUE)) |>
        ungroup()
      message("computed ", av, " (national weekly mean removed)")
    }
    if (mv %in% names(d)) {
      d[[mv]] <- scale1(d[[mv]]); fixed <- c(fixed, mv)
      message("including mobility: ", mv, "  (lag ", mob_lag, " weeks)")
    } else message("mobility requested but ", mv, " is not available; skipped")
  }

  ## ---- SEASONAL HARMONICS ---------------------------------------------------
  ## sin/cos of the week of year. REDUNDANT whenever a flexible temporal term is
  ## present: with alpha_t (RW1) and phi_it in the model the harmonics were
  ## unidentifiable in the NIMBLE fit (R-hat above 5). They are worth including
  ## only with interaction = "none" and no trend, where they become the sole
  ## representation of seasonality.
  if (harmonics) {
    hv <- intersect(c("sin1", "cos1", "sin2", "cos2"), names(d))
    if (length(hv)) {
      fixed <- c(fixed, hv)
      message("including harmonics: ", paste(hv, collapse = ", "))
      if (interaction != "none")
        message("  NOTE: harmonics alongside interaction = '", interaction,
                "' are typically unidentifiable; check the fixed effects.")
    } else message("harmonics requested but sin1/cos1 are not available")
  }

  ## Any further covariate the caller has already put in the data, for terms
  ## built outside this function -- the Florida seeding pressure, say. They are
  ## standardised like the rest so their coefficients are comparable, and a
  ## missing column fails here with a clear message rather than inside INLA's
  ## formula parser.
  if (length(extra_fixed)) {
    miss <- setdiff(extra_fixed, names(d))
    if (length(miss))
      stop("extra_fixed column(s) not in the data: ",
           paste(miss, collapse = ", "),
           "\n  add them to dat$data before fitting, e.g. with ",
           "add_florida_pressure().", call. = FALSE)
    for (v in extra_fixed) {
      n_na <- sum(!is.finite(d[[v]]))
      if (n_na) {
        ## INLA drops rows with a missing covariate, which silently changes the
        ## sample and makes WAIC incomparable with a model that lacks the term
        message("  ", v, ": ", n_na, " non-finite values set to 0 after ",
                "standardising, so the row count stays comparable")
      }
      d[[v]] <- scale1(d[[v]])
      d[[v]][!is.finite(d[[v]])] <- 0
    }
    fixed <- c(fixed, extra_fixed)
    message("including extra fixed effect(s): ",
            paste(extra_fixed, collapse = ", "))
  }

  fixed <- c(fixed, zvars)

  ## ---- random-slope indices -------------------------------------------------
  ## One index column per climate covariate, so each gets its own gamma_i.
  slope_terms <- character(0)
  if (slope_model != "none" && length(zvars)) {
    gfile_s <- write_inla_graph(W)
    for (k in seq_along(zvars)) {
      idx <- paste0("i_sl", k)
      d[[idx]] <- d$i
      slope_terms <- c(slope_terms,
        if (slope_model == "iid")
          sprintf("f(%s, %s, model = 'iid')", idx, zvars[k])
        else
          sprintf("f(%s, %s, model = 'besag', graph = '%s', scale.model = TRUE, constr = TRUE)",
                  idx, zvars[k], gfile_s))
    }
    message("state-specific climate slopes (gamma_i): ", slope_model,
            " for ", paste(zvars, collapse = ", "))
  }

  gfile <- write_inla_graph(W)

  ## ---- formula --------------------------------------------------------------
  ## The shared national trend alpha_t is optional. Dropping it is coherent ONLY
  ## if something else carries the temporal signal: the type II interaction can,
  ## since it has one term per state-week and its state-wise mean plays the role
  ## of a national curve. With a CONSTRAINED interaction (harmonic, none) there
  ## is then no shared temporal structure at all and each state is modelled
  ## independently in time -- a real change, not a simplification.
  rhs <- c("1", fixed,
           sprintf("f(i, model = 'bym2', graph = '%s', scale.model = TRUE, constr = TRUE)",
                   gfile),
           if (add_trend) "f(t, model = 'rw1', scale.model = TRUE, constr = TRUE)",
           slope_terms)
  if (!add_trend) {
    message("alpha_t OMITTED",
            if (interaction %in% c("none", "harmonic"))
              " -- and the interaction is constrained, so there is NO shared temporal term"
            else " -- the interaction now carries the national signal")
  }

  ## ---- the interaction, with control over its degrees of freedom ----------
  ## type2/type4 place one term on every (state, week) cell: S*T = 21,522 here,
  ## i.e. one per observation. The RW1 prior smooths them, but the model is
  ## close to saturated and will track almost any series. Check the EFFECTIVE
  ## number of parameters with res$fit$dic$p.eff.
  ##
  ## The alternatives below keep state-specific timing with far fewer degrees of
  ## freedom, so the fit becomes a test of structure rather than interpolation.
  if (interaction == "type2") {
    rhs <- c(rhs, "f(t2, model = 'rw1', scale.model = TRUE, constr = TRUE,
                     group = i2, control.group = list(model = 'iid'))")

  } else if (interaction == "type4") {
    rhs <- c(rhs, sprintf("f(t2, model = 'rw1', scale.model = TRUE, constr = TRUE,
                     group = i2, control.group = list(model = 'besag', graph = '%s'))",
                     gfile))

  } else if (interaction == "harmonic") {
    ## STATE-SPECIFIC SEASONAL SHAPE, 2 (or 4) parameters per state.
    ##   A_i sin(2*pi*w/52 + psi_i)  =  a_i sin(.) + b_i cos(.)
    ## so random slopes on sin and cos give each state its own seasonal
    ## AMPLITUDE and PHASE -- the travelling wave -- in 2S terms rather than S*T.
    need <- c("sin1", "cos1")
    if (!all(need %in% names(d)))
      stop("harmonic interaction needs sin1/cos1 from stcar_prepare()", call. = FALSE)
    d$i_s1 <- d$i; d$i_c1 <- d$i
    rhs <- c(rhs, "f(i_s1, sin1, model = 'iid')",
                  "f(i_c1, cos1, model = 'iid')")
    if (all(c("sin2", "cos2") %in% names(d))) {
      d$i_s2 <- d$i; d$i_c2 <- d$i
      rhs <- c(rhs, "f(i_s2, sin2, model = 'iid')",
                    "f(i_c2, cos2, model = 'iid')")
    }
    message("interaction: HARMONIC -- state-specific amplitude and phase (",
            length(states) * 4, " terms, vs ", length(states) *
            length(unique(d$t)), " for type II)")

  } else if (interaction == "harmonic_ar1") {
    ## THE FORECASTING STRUCTURE. Two terms doing different jobs:
    ##
    ##   harmonic  state-specific sin/cos of week-of-year. Deterministic in
    ##             time, so the seasonal SHAPE projects past the origin -- each
    ##             state keeps its own amplitude and phase into the forecast.
    ##
    ##   AR1       a mean-reverting local deviation. Where a random walk carries
    ##             the last value forward indefinitely, an AR1 decays toward
    ##             zero at rate rho, so a state currently running above its
    ##             seasonal average is forecast to stay above it briefly and
    ##             then return. That is what a "current conditions" term should
    ##             do; a random walk instead flattens onto the last level and
    ##             the forecast becomes persistence shifted by the horizon.
    ##
    ## Together: the climatological season carries the level, the AR1 carries
    ## recent anomaly, and both extrapolate sensibly. Pair with
    ## add_trend = FALSE -- an RW1 national trend also flattens and would drag
    ## the forecast back toward the last observed level.
    need <- c("sin1", "cos1")
    if (!all(need %in% names(d)))
      stop("harmonic_ar1 needs sin1/cos1 from stcar_prepare()", call. = FALSE)
    d$i_s1 <- d$i; d$i_c1 <- d$i
    rhs <- c(rhs, "f(i_s1, sin1, model = 'iid')",
                  "f(i_c1, cos1, model = 'iid')")
    if (all(c("sin2", "cos2") %in% names(d))) {
      d$i_s2 <- d$i; d$i_c2 <- d$i
      rhs <- c(rhs, "f(i_s2, sin2, model = 'iid')",
                    "f(i_c2, cos2, model = 'iid')")
    }
    rhs <- c(rhs, "f(t2, model = 'ar1',
                     group = i2, control.group = list(model = 'iid'))")
    message("interaction: HARMONIC + AR1 -- seasonal shape projects forward, ",
            "local deviation decays toward it")

  } else if (interaction == "harmonic_season") {
    ## WHY THE PEAKS ARE UNDERSHOT, AND WHAT THIS CHANGES
    ##
    ## `harmonic_ar1` gives each state ONE seasonal amplitude, estimated across
    ## every season it has. Forecasting a season that is running large therefore
    ## pulls toward the state's average season, and the AR1 -- which reverts to
    ## zero -- cannot hold the difference for long. That is the mechanism behind
    ## predictions that track the epidemic's timing but fall short of its peak.
    ##
    ## Here the harmonic coefficients vary by STATE AND SEASON:
    ##
    ##     a_{i,s} ~ N(a_i, sigma_a^2)
    ##
    ## so a season's amplitude is learned from its own early weeks and then
    ## projected across the rest of it, while being shrunk toward that state's
    ## usual season. Early in a season the shrinkage dominates, as it should
    ## when little has been seen; as weeks accumulate the season's own level
    ## takes over.
    ##
    ## The AR1 stays, carrying week-to-week deviation from that season's curve.
    need <- c("sin1", "cos1")
    if (!all(need %in% names(d)))
      stop("harmonic_season needs sin1/cos1 from stcar_prepare()", call. = FALSE)
    if (!"season" %in% names(d)) d$season <- season_of(d$week)

    ## one index per state-season
    d$is_idx <- as.integer(factor(paste(d$state, d$season)))
    d$is_s1 <- d$is_idx; d$is_c1 <- d$is_idx
    ## and the state-level term the season-level one shrinks toward
    d$i_s1 <- d$i; d$i_c1 <- d$i

    rhs <- c(rhs,
             "f(i_s1, sin1, model = 'iid')",
             "f(i_c1, cos1, model = 'iid')",
             "f(is_s1, sin1, model = 'iid')",
             "f(is_c1, cos1, model = 'iid')",
             "f(t2, model = 'ar1',
                group = i2, control.group = list(model = 'iid'))")
    message("interaction: HARMONIC (state and state-season) + AR1 -- ",
            "amplitude adapts to the season in progress  [",
            max(d$is_idx), " state-seasons]")

  } else if (interaction == "monthly") {
    ## Same structure as type II but on MONTHS, cutting the time dimension by
    ## roughly four while keeping each state its own trajectory.
    d$t_m  <- as.integer(factor(format(d$week, "%Y-%m")))
    d$i_m  <- d$i
    rhs <- c(rhs, "f(t_m, model = 'rw1', scale.model = TRUE, constr = TRUE,
                     group = i_m, control.group = list(model = 'iid'))")
    message("interaction: MONTHLY type II (", length(states) * max(d$t_m),
            " terms, vs ", length(states) * length(unique(d$t)), ")")

  } else if (interaction == "ar1") {
    ## AR1 is stationary and mean-reverting, so it uses far fewer EFFECTIVE
    ## degrees of freedom than a random walk, and does not wander. It also
    ## extrapolates sensibly, which matters if this is ever used to forecast.
    rhs <- c(rhs, "f(t2, model = 'ar1',
                     group = i2, control.group = list(model = 'iid'))")
    message("interaction: AR1 in time (mean-reverting; fewer effective df than RW1)")
  }

  ## ---- priors borrowed from another fit ----------------------------------
  ## Each entry is keyed by the hyperparameter label INLA reports ("Precision
  ## for i", and so on) and is spliced into the matching f() term. `rhs` is a
  ## character vector at this point, so the prior is injected textually.
  ##
  ## A label that matches no term in THIS model is reported rather than
  ## ignored: a run where nothing matched produces exactly the same fit as one
  ## with default priors, and without the message there is no way to tell the
  ## two apart. That is worth guarding, because the labels depend on the
  ## interaction -- a harmonic_season source fit carries is_s1 and is_c1,
  ## which a harmonic_ar1 model does not have.
  if (!is.null(hyper_priors)) {
    used <- character(0)
    for (nm in names(hyper_priors)) {
      v <- sub("^Precision for ", "", nm)
      v <- sub("[ (].*$", "", v)
      ## rhs entries may carry leading whitespace and newlines
      j <- which(grepl(paste0("^\\s*f\\(\\s*", v, "\\s*,"), rhs))
      if (!length(j)) next
      ## if this term already carries a hyper list from an earlier entry, the
      ## new one has to be merged into it rather than appended as a second
      ## argument, which INLA would reject
      if (grepl("hyper\\s*=\\s*list\\(", rhs[j])) {
        add <- character(0)
        hp2 <- hyper_priors[[nm]]
        if (!is.null(hp2$rho))
          add <- sprintf("rho = list(prior = '%s', param = c(%g, %g))",
                         hp2$rho$prior, hp2$rho$param[1], hp2$rho$param[2])
        if (!is.null(hp2$prec) && !grepl("prec\\s*=", rhs[j]))
          add <- c(add, sprintf("prec = list(prior = '%s', param = c(%g, %g))",
                                hp2$prec$prior, hp2$prec$param[1],
                                hp2$prec$param[2]))
        if (length(add)) {
          rhs[j] <- sub("hyper\\s*=\\s*list\\(",
                        paste0("hyper = list(", paste(add, collapse = ", "),
                               ", "), rhs[j])
          used <- c(used, nm)
        }
        next
      }
      hp <- hyper_priors[[nm]]
      ## a prior may be on the precision OR on the AR(1) correlation; both go
      ## inside the same `hyper` list, keyed by the name INLA expects
      parts <- character(0)
      if (!is.null(hp$prec))
        parts <- c(parts, sprintf("prec = list(prior = '%s', param = c(%g, %g))",
                                  hp$prec$prior, hp$prec$param[1],
                                  hp$prec$param[2]))
      if (!is.null(hp$rho))
        parts <- c(parts, sprintf("rho = list(prior = '%s', param = c(%g, %g))",
                                  hp$rho$prior, hp$rho$param[1],
                                  hp$rho$param[2]))
      if (!length(parts)) next
      txt <- paste0("hyper = list(", paste(parts, collapse = ", "), ")")
      ## close on the LAST bracket of the term, not the first
      rhs[j] <- sub("\\)\\s*$", paste0(", ", txt, ")"), rhs[j])
      used <- c(used, nm)
    }
    message("borrowed priors applied to: ",
            if (length(used)) paste(used, collapse = ", ") else "NOTHING")
    un <- setdiff(names(hyper_priors), used)
    if (length(un))
      message("  no matching term, so NOT applied: ", paste(un, collapse = ", "))
    if (!length(used))
      warning("no borrowed prior matched a term in this model; the fit is ",
              "identical to one with default priors", call. = FALSE)
  }

  form <- stats::as.formula(paste("y_inla ~", paste(rhs, collapse = " + ")))
  message("\nformula:\n  ", deparse(form, width.cutoff = 200))
  message("family: ", fam)

  t0 <- Sys.time()
  fit <- INLA::inla(form, family = fam, data = d, E = d$encounters,
                    control.family = cf,
                    ## INLA reports only the quantile levels asked for here.
                    ## The default three are enough for a 95% interval; a
                    ## forecast hub wanting 23 must have them requested at fit
                    ## time, because a summary cannot be un-summarised later.
                    quantiles = if (is.null(quantiles))
                      c(0.025, 0.5, 0.975) else sort(unique(quantiles)),
                    control.predictor = list(compute = TRUE, link = 1),
                    ## `light = TRUE` drops everything a nowcast does not use.
                    ## config = TRUE stores the whole GMRF configuration so that
                    ## inla.posterior.sample() can run, and it dominates the size
                    ## of the fitted object; dic and waic add more. A rolling
                    ## nowcast needs only the fitted values, so carrying them
                    ## makes each fit several times larger for no benefit -- the
                    ## usual cause of memory exhaustion partway through a run.
                    control.compute = if (light)
                      list(dic = FALSE, waic = FALSE, config = FALSE)
                    else list(dic = TRUE, waic = TRUE, config = TRUE),
                    verbose = verbose)
  message(sprintf("INLA finished in %.1f min",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  d$fit_mean <- fit$summary.fitted.values[, "mean"]
  d$mu_med <- fit$summary.fitted.values[, "0.5quant"]  * d$encounters
  d$mu_lo  <- fit$summary.fitted.values[, "0.025quant"] * d$encounters
  d$mu_hi  <- fit$summary.fitted.values[, "0.975quant"] * d$encounters

  ## EXTRA QUANTILES, for a forecast hub that wants more than three.
  ## They have to be asked for at fit time -- INLA reports only the levels
  ## given in `quantiles`, and there is no way to recover the rest afterwards
  ## from a summary. Interpolating 23 levels through three knots would be
  ## wrong exactly in the tails, which is where a proper score is most
  ## sensitive, so the levels are requested rather than reconstructed.
  if (!is.null(quantiles)) {
    qn <- sprintf("%squant", quantiles)
    have <- intersect(qn, colnames(fit$summary.fitted.values))
    for (q in have)
      d[[paste0("q", sub("quant$", "", q))]] <-
        fit$summary.fitted.values[, q] * d$encounters
    miss <- setdiff(qn, have)
    if (length(miss))
      warning(length(miss), " requested quantile(s) not returned by INLA: ",
              paste(head(miss, 5), collapse = ", "), call. = FALSE)
  }

  if (!is.null(fit$dic$p.eff))
    message(sprintf("effective number of parameters (p.eff) = %.0f  of %d observations",
                    fit$dic$p.eff, nrow(d)))

  if (light) {
    ## keep only what downstream code reads; the rest of an INLA object is
    ## large and is never touched by the nowcast
    fit <- list(summary.fitted.values = fit$summary.fitted.values,
                summary.fixed  = fit$summary.fixed,
                summary.hyperpar = fit$summary.hyperpar,
                summary.random = fit$summary.random[
                  intersect(names(fit$summary.random), c("i", "t"))])
  }

  out <- list(fit = fit, data = d, states = states, formula = form,
              family = fam, interaction = interaction, add_trend = add_trend,
              cen_lo = cen_lo, cen_hi = cen_hi,
              slope_model = slope_model, zvars = zvars, censored = censored)
  class(out) <- c("stcar_inla3", "list")
  out
}


## ---------------------------------------------------------------------------
## diagnostics for the new pieces
## ---------------------------------------------------------------------------
#' Is the NB dispersion identified, or has the interaction absorbed it?
#' Is the AR1 actually mean-reverting?
#'
#' An AR1 with rho near 1 IS a random walk: it carries the last deviation
#' forward instead of decaying back to the seasonal curve, and the forecast
#' becomes persistence again. That would explain a prediction correlating 0.99
#' with the value at the origin despite the seasonal terms being present. Values
#' around 0.7 to 0.9 mean the deviation genuinely decays.
check_rho <- function(res) {
  h <- res$fit$summary.hyperpar
  i <- grep("^Rho", rownames(h))
  if (!length(i)) { cat("\nno AR1 term in this fit\n"); return(invisible(NULL)) }
  q <- h[i[1], c("0.025quant", "0.5quant", "0.975quant")]
  cat(sprintf("\nAR1 rho = %.3f  [%.3f, %.3f]\n", q[[2]], q[[1]], q[[3]]))
  ## the deviation decays by rho each week, so this is how long it persists
  hl <- log(0.5) / log(max(min(q[[2]], 0.999), 1e-6))
  cat(sprintf("half-life of a deviation: %.1f weeks\n", hl))
  if (q[[2]] > 0.95)
    cat("  rho above 0.95: the AR1 is behaving like a random walk, so the\n",
        "  forecast will stay close to the last observed value however good\n",
        "  the seasonal terms are.\n", sep = "")
  else
    cat("  the deviation decays back toward the seasonal curve, which is what\n",
        "  lets the forecast differ from persistence at longer horizons.\n",
        sep = "")
  invisible(q)
}


check_theta <- function(res) {
  h <- res$fit$summary.hyperpar
  th <- grep("size|dispersion", rownames(h), value = TRUE, ignore.case = TRUE)
  if (!length(th)) { cat("\nno NB dispersion in this fit (Poisson family)\n"); return(invisible(NULL)) }
  q <- h[th[1], c("0.025quant", "0.5quant", "0.975quant")]
  cat(sprintf("\nNB dispersion: %.2f  [%.2f, %.2f]\n", q[[2]], q[[1]], q[[3]]))
  if (q[[2]] > 100 || q[[3]] / max(q[[1]], 1e-9) > 20)
    cat("!! very large or very wide: the interaction phi_it is absorbing the\n",
        "   overdispersion and theta is not identified. Use cenpoisson, or drop\n",
        "   the interaction -- not both.\n", sep = "")
  else
    cat("   plausible: the NB dispersion appears identified alongside phi_it.\n")
  invisible(q)
}


#' State-specific climate slopes gamma_i.
#' Returns the common slope beta, the state deviations, and beta + gamma_i.
climate_slopes <- function(res) {
  if (!length(res$zvars) || res$slope_model == "none") {
    cat("\nno state-specific climate slopes in this fit\n"); return(invisible(NULL))
  }
  out <- list()
  for (k in seq_along(res$zvars)) {
    v   <- res$zvars[k]
    idx <- paste0("i_sl", k)
    r   <- res$fit$summary.random[[idx]]
    if (is.null(r)) next
    b   <- res$fit$summary.fixed[v, "0.5quant"]
    tab <- data.frame(state = res$states,
                      common_slope = b,
                      gamma_i = r[["0.5quant"]],
                      slope_i = b + r[["0.5quant"]],
                      lo = b + r[["0.025quant"]],
                      hi = b + r[["0.975quant"]],
                      row.names = NULL)
    tab$differs <- ifelse(r[["0.025quant"]] > 0 | r[["0.975quant"]] < 0, "yes", "no")
    out[[v]] <- tab
    cat("\n--- ", v, " ---\n", sep = "")
    cat(sprintf("common slope beta = %.3f\n", b))
    cat(sprintf("state slopes range %.3f to %.3f; %d of %d differ credibly\n",
                min(tab$slope_i), max(tab$slope_i), sum(tab$differs == "yes"),
                nrow(tab)))
    print(utils::head(tab[order(tab$slope_i), c("state","slope_i","lo","hi","differs")], 5))
  }
  invisible(out)
}


#' Map the state-specific slopes: where does climate matter most?
plot_climate_slopes <- function(res, var = NULL,
                                file = NULL) {
  suppressPackageStartupMessages(library(ggplot2))
  sl <- climate_slopes(res)
  if (is.null(sl) || !length(sl)) return(invisible(NULL))
  var <- var %||% names(sl)[1]
  tab <- sl[[var]]

  p <- ggplot(tab, aes(reorder(state, slope_i), slope_i)) +
    geom_col(aes(fill = differs)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .3, linewidth = .3) +
    geom_hline(yintercept = tab$common_slope[1], linetype = "dashed") +
    coord_flip() +
    scale_fill_manual(values = c(no = "grey70", yes = "#b2182b"),
                      name = "differs from\ncommon slope") +
    labs(title = paste0("State-specific slope for ", var),
         subtitle = "dashed line = common slope; bars = beta + gamma_i (log scale per 1 SD)",
         x = NULL, y = "slope") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 6))
  if (!is.null(file)) { ggsave(file, p, width = 9, height = 11, dpi = 150)
                        message("wrote ", file) }
  p
}



## ---------------------------------------------------------------------------
## compare interaction structures
## ---------------------------------------------------------------------------
#' Fit the same model with several interaction structures and compare.
#'
#' The fitted-curve figure cannot distinguish these -- a saturated model tracks
#' the data whatever else is in it -- so the comparison has to be made on
#' information criteria and effective degrees of freedom.
#'
#' Read the result as: if `none` is far worse than the others, the interaction
#' matters. If `harmonic` is close to `type2` at 1/200th the parameters, the
#' extra flexibility of type II is buying interpolation rather than structure.
#'
#' @param which interaction structures to compare
compare_interactions <- function(dat,
                                 which = c("none", "harmonic", "monthly", "type2"),
                                 climate = NULL, use_pop = FALSE,
                                 family = "cenpoisson", subset_states = NULL,
                                 verbose = FALSE) {
  rows <- list(); fits <- list()
  for (w in which) {
    message("\n===== interaction = ", w, " =====")
    r <- try(stcar_fit_inla3(dat, climate = climate, use_pop = use_pop,
                             slope_model = "none", family = family,
                             interaction = w, censored = TRUE,
                             subset_states = subset_states, verbose = verbose),
             silent = TRUE)
    if (inherits(r, "try-error")) {
      warning("interaction = ", w, " failed: ",
              conditionMessage(attr(r, "condition")), call. = FALSE)
      next
    }
    fits[[w]] <- r
    rows[[w]] <- data.frame(
      interaction = w,
      DIC   = r$fit$dic$dic,
      WAIC  = r$fit$waic$waic,
      p_eff = r$fit$dic$p.eff,
      n_obs = nrow(r$data))
  }
  tab <- do.call(rbind, rows)
  if (!is.null(tab)) {
    tab <- tab[order(tab$WAIC), ]
    tab$dWAIC <- round(tab$WAIC - min(tab$WAIC), 1)
    tab$p_eff_per_obs <- round(tab$p_eff / tab$n_obs, 3)
    cat("\n=========== INTERACTION COMPARISON ===========\n")
    print(tab, row.names = FALSE)
    cat("\np_eff_per_obs near 1 means one effective parameter per observation,\n",
        "i.e. the model is interpolating rather than describing structure.\n", sep = "")
  }
  invisible(list(table = tab, fits = fits))
}
