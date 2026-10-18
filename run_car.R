## run_car.R — fit and backtest the CAR model on the NHSN RSV panel
##   source("run_car.R")      # run from the folder containing these files

suppressPackageStartupMessages({ library(dplyr); library(readr) })
for (f in c("helpers.R", "stcar_covariates.R", "stcar_prep.R", "stcar_inla3.R",
            "nowcast_inla.R", "hub_data.R", "florida_seed.R")) source(f)

## data: NHSN RSV admissions from the CDC RSV Forecast Hub (cached after the
## first download); population as the offset
hub <- hub_target_data("hosp", cache = "data_hub", verbose = FALSE)
hp  <- hub_to_panel(hub, drop_locations = c("US", "72"), verbose = FALSE)
st  <- sort(unique(hp$data$state))
pop <- hub$locations |> filter(location_name %in% st) |>
  { \(x) setNames(as.numeric(x$population), x$location_name) }()
dat <- stcar_panel(hp$data |> mutate(y_raw = as.integer(y_raw)),
                   states = st, pop = pop, verbose = FALSE)
dat <- add_florida_pressure(dat, seed = "Florida", weight = "exp", lag = 4,
                            rho_km = 1500, verbose = FALSE)

## one fit on all the data
fit <- stcar_fit_inla3(dat, interaction = "harmonic_ar1", add_trend = FALSE,
                       family = "poisson", censored = FALSE,
                       extra_fixed = "fl_exp_lag4",
                       climate = NULL, use_pop = FALSE, mobility = "none",
                       slope_model = "none")
print(fit$fit$summary.hyperpar)

## rolling backtest over the 2025-26 season, forecasts 1-4 weeks ahead
fc <- rolling_nowcast(dat, valid_start = as.Date("2025-07-01"),
                      valid_end = as.Date("2026-06-27"), horizons = 1:4,
                      refit_every = 2, interaction = "harmonic_ar1",
                      add_trend = FALSE, family = "poisson", censored = FALSE,
                      extra_fixed = "fl_exp_lag4",
                      climate = NULL, use_pop = FALSE, mobility = "none",
                      slope_model = "none", checkpoint_dir = "car_origins")
write_csv(fc, "car_forecasts.csv")

## score: weighted interval score at the 95% interval, by horizon
fc |> left_join(dat$data |> select(state, week, truth = y_raw),
                by = c("state", "week")) |>
  filter(is.finite(truth)) |>
  mutate(WIS = .wis1(truth, mu_med, mu_lo, mu_hi, 0.05)) |>
  group_by(horizon) |>
  summarise(WIS = mean(WIS), coverage = mean(truth >= mu_lo & truth <= mu_hi)) |>
  print()
