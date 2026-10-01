
# Load libraries and functions --------------------------------------------

source("scripts/0. Setup/0.0 Libraries and General Functions.R")

# Set default backend to cmdstanr and reduce verbosity
options(brms.backend = "cmdstanr")
Sys.setenv(CMDSTANR_NO_VERBOSITY = 1)

# Generate simulation design matrices ----------------------------------------------------------------

## Setup Simulations 1 to 3: Basic Concept -------------

# Setup the grid for Simulations 1 through 3
simulation1 = scenario_grid(20, c(.80, .90), c(0.15, 0.3), 0.2, c(.015, .06, .12), c(0.5, 1.0), -.5, 0.0, 200)
simulation2 = scenario_grid(c(10, 20, 40), .85, 0.2, c(0.1, 0.3), .03, 1.0, -.5, 0.0, 200)
simulation3 = scenario_grid(20, .85, 0.2, 0.2, .03, 1.0, -.5, c(0.0, 0.1, 0.3), c(200, 300))

## Setup Simulation 4: Systematic Study-level Prior Misspecification -------------

# Setup the grid for Simulations 4a and 4b; note that only 4a is reported in-text
simulation4a = scenario_grid_sim4(
  sn_levels     = 20,
  accuracy      = 0.85,
  sigma_levels  = 0.2,
  omega_levels  = 0.2,
  mu_prev_vals  = c(0.015, 0.06, 0.12),
  tau           = 1.0,
  rho           = -0.5,
  sigma_bias    = 0.3,
  kappa_prior   = 200,
  delta_se_vals = c(-0.10, 0, 0.10),
  delta_sp_vals = c(-0.10, 0, 0.10)
)

simulation4b = scenario_grid_sim4(
  sn_levels     = 20,
  accuracy      = 0.85,
  sigma_levels  = 0.2,
  omega_levels  = 0.2,
  mu_prev_vals  = c(0.015, 0.06, 0.12),
  tau           = 1.0,
  rho           = -0.5,
  sigma_bias    = 0.3,
  kappa_prior   = 50,
  delta_se_vals = c(-0.10, 0, 0.10),
  delta_sp_vals = c(-0.10, 0, 0.10)
)

## Setup Simulation 5: Common Prior Center -------------

# Setup the grid for Simulations 5a and 5b; note that only 5a is reported in-text
simulation5a = scenario_grid_sim5(
  sn_levels    = 20,
  accuracy     = c(0.80, 0.90),
  sigma_levels = 0.2,
  omega_levels = 0.2,
  mu_prev_vals = c(0.015, 0.06, 0.12),
  tau          = 1.0,
  rho          = -0.5,
  sigma_bias   = 0,
  kappa_prior  = 200
)

simulation5b = scenario_grid_sim5(
  sn_levels    = 20,
  accuracy     = c(0.80, 0.90),
  sigma_levels = 0.2,
  omega_levels = 0.2,
  mu_prev_vals = c(0.015, 0.06, 0.12),
  tau          = 1.0,
  rho          = -0.5,
  sigma_bias   = 0.3,
  kappa_prior  = 200
)


## Setup Simulation 6: Gold-Standard Anchoring -----------------------------

simulation6 = scenario_grid_sim6(
  sn_levels      = 20,
  accuracy       = 0.85,
  sigma_levels   = 0.2,
  omega_levels   = 0.2,
  mu_prev_vals   = c(0.015, 0.06, 0.12),
  tau            = 1.0,
  rho            = -0.5,
  sigma_bias     = 0.3,
  kappa_prior    = 200,
  delta_sp_vals  = c(-0.10, 0, 0.10),
  prop_gold_vals = c(0.10, 0.25, 0.50)
)

## Setup Simulation 7: Moderator Detection ---------------------------------

simulation7 = scenario_grid_sim7(
  sn_levels      = 40,
  accuracy       = 0.85,
  sp_offset      = 0.10,
  sigma_levels   = 0.2,
  omega_levels   = 0.2,
  mu_prev_vals   = 0.03,
  tau            = 1.0,
  rho            = -0.5,
  kappa_prior    = 100,
  slope_vals     = c(0, 0.5, 1.0),
  prop_gold_vals = c(0.10, 0.25, 0.50),
  prior_se_vals  = c(0.85, 0.85),
  prior_sp_vals  = c(0.75, 0.85)
)


## Setup Simulation 8: Se/Sp Estimation ------------------------------------

simulation8 = scenario_grid_sim8(
  sn_levels      = c(20, 40),
  mean_se_vals   = 0.85,
  mean_sp_vals   = c(0.70, 0.80, 0.90),
  sigma_levels   = 0.2,
  omega_levels   = 0.2,
  mu_prev_vals   = c(0.015, 0.06, 0.12),
  tau            = 1.0,
  rho            = -0.5,
  sigma_bias     = 0,
  kappa_prior    = 100,
  prop_gold_vals = c(0.10, 0.25, 0.50)
)

# Generate or Load Sim Data -------------------------------------------------------

#
# UNCOMMENT TO GENERATE AND SAVE THE DATA
#

# sim1_sims = simulate_all_scenarios(R = 200, grid = simulation1)
# sim2_sims = simulate_all_scenarios(R = 200, grid = simulation2)
# sim3_sims = simulate_all_scenarios(R = 200, grid = simulation3)
# sim4a_sims = simulate_all_scenarios_sim4(R = 200, grid = simulation4a)
# sim4b_sims = simulate_all_scenarios_sim4(R = 200, grid = simulation4b)
# sim5a_sims = simulate_all_scenarios_sim5(R = 200, grid = simulation5a)
# sim5b_sims = simulate_all_scenarios_sim5(R = 200, grid = simulation5b)
# sim6_sims = simulate_all_scenarios_sim6(R = 200, grid = simulation6)
# sim7_sims = simulate_all_scenarios_sim7(R = 200, grid = simulation7)
# sim8_sims = simulate_all_scenarios_sim8(R = 200, grid = simulation8)
#
# saveRDS(sim1_sims, 'simulation_data/sim1/sim1_dat.rds')
# saveRDS(sim2_sims, 'simulation_data/sim2/sim2_dat.rds')
# saveRDS(sim3_sims, 'simulation_data/sim3/sim3_dat.rds')
# saveRDS(sim4a_sims, 'simulation_data/sim4a/sim4a_dat.rds')
# saveRDS(sim4b_sims, 'simulation_data/sim4b/sim4b_dat.rds')
# saveRDS(sim5a_sims, 'simulation_data/sim5a/sim5a_dat.rds')
# saveRDS(sim5b_sims, 'simulation_data/sim5b/sim5b_dat.rds')
# saveRDS(sim6_sims, 'simulation_data/sim6/sim6_dat.rds')
# saveRDS(sim7_sims, 'simulation_data/sim7/sim7_dat.rds')
# saveRDS(sim8_sims, 'simulation_data/sim8/sim8_dat.rds')

sim1_sims = readRDS('simulation_data/sim1/sim1_dat.rds')
sim2_sims = readRDS('simulation_data/sim2/sim2_dat.rds')
sim3_sims = readRDS('simulation_data/sim3/sim3_dat.rds')
sim4a_sims = readRDS('simulation_data/sim4a/sim4a_dat.rds')
sim4b_sims = readRDS('simulation_data/sim4b/sim4b_dat.rds')
sim5a_sims = readRDS('simulation_data/sim5a/sim5a_dat.rds')
sim5b_sims = readRDS('simulation_data/sim5b/sim5b_dat.rds')
sim6_sims = readRDS('simulation_data/sim6/sim6_dat.rds')
sim7_sims = readRDS('simulation_data/sim7/sim7_dat.rds')
sim8_sims = readRDS('simulation_data/sim8/sim8_dat.rds')

# Note: Individual replicates can be accessed via dat_raw = sim1_sims[[1]]$replicates[[1]]

# Run all simulations -----------------------------------------------------

#
# UNCOMMENT ONLY IF YOU INTEND TO RUN THESE: They are very time consuming
# and have been commented out to avoid accidentally running this code. Once run,
# the data are saved and can be accessed via the following script.
#


## Run Simulation 1 --------------------------------------------------------
#
# run_all_rds(sim1_sims, iter = 6000, adapt_delta = 0.9999,
#                 dir='simulation_data/sim1/', stan_m7 = stan_m7_unbounded,
#                 is_bounded = FALSE)
#


## Run Simulation 2 --------------------------------------------------------
#
# run_all_rds(sim2_sims, iter = 6000, adapt_delta = 0.9999,
#                 dir='simulation_data/sim2/', stan_m7 = stan_m7_unbounded,
#                 is_bounded = FALSE)
#


## Run Simulation 3 --------------------------------------------------------
#
# run_all_rds(sim3_sims, iter = 6000, adapt_delta = 0.9999,
#                 dir='simulation_data/sim3/', stan_m7 = stan_m7_unbounded,
#                 is_bounded = FALSE)
#


## Run Simulation 4 --------------------------------------------------------
#
# run_all_rds_m1m6(sim4a_sims, dir = 'simulation_data/sim4a/',
#                  fit_fn = fit_m1_m6_sim4,
#                  iter = 6000, adapt_delta = 0.9999)
#
# run_all_rds_m1m6(sim4b_sims, dir = 'simulation_data/sim4b/',
#                  fit_fn = fit_m1_m6_sim4,
#                  iter = 6000, adapt_delta = 0.9999)
#


## Run Simulation 5 --------------------------------------------------------
#
# run_all_rds_m1m6(sim5a_sims, dir = 'simulation_data/sim5a/',
#                  fit_fn = fit_m1_m6_sim5,
#                  iter = 6000, adapt_delta = 0.9999)
#
# run_all_rds_m1m6(sim5b_sims, dir = 'simulation_data/sim5b/',
#                  fit_fn = fit_m1_m6_sim5,
#                  iter = 6000, adapt_delta = 0.9999)

## Run Simulation 6 --------------------------------------------------------
#
# run_all_rds_m1m6(sim6_sims, dir = 'simulation_data/sim6/',
#                  fit_fn = fit_m1_m6_sim6,
#                  iter = 6000, adapt_delta = 0.9999)

## Run Simulation 7 --------------------------------------------------------
#
# run_all_rds_sim7(sim7_sims, dir = 'simulation_data/sim7/',
#                  iter = 6000, adapt_delta = 0.9999)

## Run Simulation 8 --------------------------------------------------------
#
# run_sim8(sim8_sims, dir = "simulation_data/sim8/",
#          iter = 6000, adapt_delta = 0.9999, seed_base = 20250321)


# Repair Simulation Runs --------------------------------------------------


## Repair Simulation 1 -----------------------------------------------------
#
# repair_cycle(sim1_sims, iter = 12000, adapt_delta = 0.9999, step_size=0.002,
#                  dir='simulation_data/sim1/', cor_dir='simulation_data/bad_runs/sim1/',
#                  cycles=3)
#


## Repair Simulation 2 -----------------------------------------------------
#
# repair_cycle(sim2_sims, iter = 12000, adapt_delta = 0.9999, step_size=0.002,
#                  dir='simulation_data/sim2/', cor_dir='simulation_data/bad_runs/sim2/',
#                  cycles=3)
#


## Repair Simulation 3 -----------------------------------------------------
#
# repair_cycle(sim3_sims, iter = 12000, adapt_delta = 0.9999, step_size=0.002,
#                  dir='simulation_data/sim3/', cor_dir='simulation_data/bad_runs/sim3/',
#                  cycles=3)
#


## Repair Simulation 4 -----------------------------------------------------
#
# repair_cycle_m1m6(sim4a_sims, fit_fn = fit_m1_m6_sim4,
#                   iter = 12000, adapt_delta = 0.9999, step_size = 0.002,
#                   dir = 'simulation_data/sim4a/',
#                   cor_dir = 'simulation_data/bad_runs/sim4a/',
#                   cycles = 3)
#
# repair_cycle_m1m6(sim4b_sims, fit_fn = fit_m1_m6_sim4,
#                   iter = 12000, adapt_delta = 0.9999, step_size = 0.002,
#                   dir = 'simulation_data/sim4b/',
#                   cor_dir = 'simulation_data/bad_runs/sim4b/',
#                   cycles = 3)
#

## Repair Simulation 5 ------------------------------------------------------------
#
# repair_cycle_m1m6(sim5a_sims, fit_fn = fit_m1_m6_sim5,
#                   iter = 12000, adapt_delta = 0.9999, step_size = 0.002,
#                   dir = 'simulation_data/sim5a/',
#                   cor_dir = 'simulation_data/bad_runs/sim5a/',
#                   cycles = 3)
#
# repair_cycle_m1m6(sim5b_sims, fit_fn = fit_m1_m6_sim5,
#                   iter = 12000, adapt_delta = 0.9999, step_size = 0.002,
#                   dir = 'simulation_data/sim5b/',
#                   cor_dir = 'simulation_data/bad_runs/sim5b/',
#                   cycles = 3)


## Repair Simulation 6 -----------------------------------------------------
#
# repair_cycle_m1m6(sim6_sims, fit_fn = fit_m1_m6_sim6,
#                   iter = 12000, adapt_delta = 0.9999, step_size = 0.002,
#                   dir = 'simulation_data/sim6/',
#                   cor_dir = 'simulation_data/bad_runs/sim6/',
#                   cycles = 3)

## Repair Simulation 7 -----------------------------------------------------
#
# repair_cycle_sim7(sim7_sims,
#                   iter = 12000, adapt_delta = 0.9999, step_size = 0.002,
#                   dir = 'simulation_data/sim7/',
#                   cor_dir = 'simulation_data/bad_runs/sim7/',
#                   cycles = 3)



## Repair Simulation 8 -----------------------------------------------------
#
# repair_cycle_sim8(sim8_sims,
#                   iter = 12000, adapt_delta = 0.9999, step_size = 0.002,
#                   dir = 'simulation_data/sim8/',
#                   cor_dir = 'simulation_data/bad_runs/sim8/', seed_base = 20250321,
#                   cycles = 3)

# Reset Parameters --------------------------------------------------------

Sys.setenv(CMDSTANR_NO_VERBOSITY = 0)   # resets Stan output
