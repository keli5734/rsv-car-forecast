# RSV CAR forecasting model

Weekly NHSN RSV hospital admissions, 51 states, forecast 1–4 weeks ahead.

## Files
- `run_car.R` — **run this**: downloads the data, fits once, runs the
  backtest and scores it. `source("run_car.R")` from this folder.
- `stcar_inla3.R` — the model: `stcar_fit_inla3()`. Only the
  `interaction = "harmonic_ar1"` branch is in use.
- `nowcast_inla.R` — forecasting: `nowcast_one()` hides weeks after the
  origin and refits; `rolling_nowcast()` runs the backtest.
- `stcar_prep.R` — `stcar_panel()` builds the indices, harmonics and offset.
- `florida_seed.R` — the Florida seeding covariate.
- `hub_data.R`, `stcar_covariates.R`, `helpers.R` — data download and helpers.

## Specification
    stcar_fit_inla3(dat, interaction = "harmonic_ar1", add_trend = FALSE,
                    family = "poisson", censored = FALSE,
                    extra_fixed = "fl_exp_lag4")

y_it ~ Poisson(mu_it)

log mu_it = log P_i + alpha + b_i
            + sum_{q=1,2} [ a_qi sin(2 pi q w_t / 52) + c_qi cos(2 pi q w_t / 52) ]
            + phi_it + beta * x_{i,t-4}

- P_i state population (offset); b_i BYM2 on a distance graph
- a_qi, c_qi state-specific harmonic coefficients, iid normal
- phi_it AR(1) within each state
- x_it Florida incidence discounted by distance, lagged 4 weeks

