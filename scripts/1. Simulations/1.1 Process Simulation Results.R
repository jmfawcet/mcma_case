
# Load libraries and functions --------------------------------------------

source("scripts/0. Setup/0.0 Libraries and General Functions.R")
# source("scripts/1. Simulations/1.0 Run Simulations.R") # Only uncomment if you want to run the simulations or are certain the run functions are commented out

# Summarize Simulation 1 ------------------------------------------------------------

sim1_sum = process_simulation(
  dir          = "simulation_data/sim1/",
  model_filter   = c("M1", "M2", "M3", "M4", "M5", "M6"),
  group_vars   = c("mean_se", "sd_se_logit", "mu_prevalence", "tau_logit"),
  exclude_models = "M7",
  exclude_bad = TRUE
)

# Summarize Simulation 2 ------------------------------------------------------------

sim2_sum = process_simulation(
  dir        = "simulation_data/sim2/",
  model_filter   = c("M1", "M2", "M3", "M4", "M5", "M6"),
  group_vars = c("sn", "omega_se"),
  exclude_models = "M7",
  exclude_bad = TRUE
)

# Simulation 3 ------------------------------------------------------------

sim3_sum = process_simulation(
  dir        = "simulation_data/sim3/",
  model_filter   = c("M1", "M2", "M3", "M4", "M5", "M6"),
  group_vars = c("sigma_bias", "kappa_prior"),
  exclude_models = "M7",
  exclude_bad = TRUE
)

# Simulation 4a ------------------------------------------------------------

sim4a_sum = process_simulation(
  dir          = "simulation_data/sim4a/",
  group_vars   = c("delta_se", "delta_sp", "mu_prevalence"),
  exclude_bad = TRUE
)

# Simulation 4b ------------------------------------------------------------

sim4b_sum = process_simulation(
  dir          = "simulation_data/sim4b/",
  group_vars   = c("delta_se", "delta_sp", "mu_prevalence"),
  exclude_bad = TRUE
)

# Simulation 5a ------------------------------------------------------------

sim5a_sum = process_simulation(
  dir          = "simulation_data/sim5a/",
  group_vars   = c("mean_se", "mu_prevalence"),
  exclude_bad = TRUE
)

# Simulation 5b ------------------------------------------------------------

sim5b_sum <- process_simulation(
  dir          = "simulation_data/sim5b/",
  group_vars   = c("mean_se", "mu_prevalence"),
  exclude_bad = TRUE
)

# Simulation 6 ------------------------------------------------------------

sim6_sum = process_simulation(
  dir          = "simulation_data/sim6/",
  group_vars   = c("prop_gold", "delta_sp", "mu_prevalence"),
  exclude_bad = TRUE
)

# Simulation 7 Analysis ----------------------------------------------------

sim7_sum <- process_simulation(
  dir = "simulation_data/sim7/",
  group_vars = c("slope_true", "prop_gold", "prior_sp"),
  model_filter = c("M1", "Gold", "KnownSeSp", "AllGold", "M6", "M8"),
  summary_type = "moderator",
  exclude_bad = TRUE
)

# Simulation 8 Analysis ---------------------------------------------------

sim8_sum = process_sim8(
  dir        = "simulation_data/sim8/",
  group_vars = c("sn","mu_prevalence", "mean_sp", "prop_gold"),
  exclude_bad = TRUE
)

sim8_sum$summary %>%
  select(sn, mu_prevalence, mean_sp, prop_gold, m6d_prev_coverage, m6j_prev_coverage, kwga_prev_bias, m6d_prev_bias, m6j_prev_bias, joint_prev_coverage, joint_prev_bias, n=n) %>%
  print(n=1000)
100*sum(sim8_sum$summary$n)/10800

sim8_sum$summary %>%
  select(sn, mu_prevalence, mean_sp, prop_gold, disc_sp_bias_wt, disc_se_bias_wt, joint_sp_bias, joint_se_bias, joint_sp_coverage, joint_se_coverage, n=n) %>%
  print(n=1000)


# Amalgamate into Excel ---------------------------------------------------

# Create a new workbook
wb = createWorkbook()

# Add each data frame as a worksheet
addWorksheet(wb, "Simulation1")
writeData(wb, "Simulation1", sim1_sum$summary)

addWorksheet(wb, "Simulation2")
writeData(wb, "Simulation2", sim2_sum$summary)

addWorksheet(wb, "Simulation3")
writeData(wb, "Simulation3", sim3_sum$summary)

addWorksheet(wb, "Simulation4a")
writeData(wb, "Simulation4a", sim4a_sum$summary)

addWorksheet(wb, "Simulation4b")
writeData(wb, "Simulation4b", sim4b_sum$summary)

addWorksheet(wb, "Simulation5a")
writeData(wb, "Simulation5a", sim5a_sum$summary)

addWorksheet(wb, "Simulation5b")
writeData(wb, "Simulation5b", sim5b_sum$summary)

addWorksheet(wb, "Simulation6")
writeData(wb, "Simulation6", sim6_sum$summary)

addWorksheet(wb, "Simulation7")
writeData(wb, "Simulation7", sim7_sum$summary)

addWorksheet(wb, "Simulation8")
writeData(wb, "Simulation8", sim8_sum$summary)

# Save the workbook
saveWorkbook(wb, file = "simulation_data/summaries/sim_summaries.xlsx", overwrite = TRUE)

