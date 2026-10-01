# ============================================================================
# Simulation 4 & 5: Systematic Prior Misspecification
# ============================================================================
#
# Simulation 4: Systematic over/underestimation of Se/Sp priors
#   - Adds a fixed shift (delta_se, delta_sp) to per-measure prior centres
#   - Also includes random perturbation (sigma_bias = 0.3)
#   - Only M1 and M6 are fit
#   - 4a uses a more concentrated kappa than 4b
#
# Simulation 5: Common prior centre for all measures
#   - Uses population-level muSe/muSp as the prior centre for EVERY measure
#   - Bias is 'variable' because each measure's true accuracy differs from
#     the common centre
#   - Only M1 and M6 are fit
#   - 5a excludes the per measure perturbation whereas 5b re-adds it (5a is more realistic)
#
# ============================================================================


# ---- Scenario Grids -------------------------------------------------------

# Scenario grid for Simulation 4
# delta_se / delta_sp: systematic shift applied to prior centres (prob scale)
scenario_grid_sim4 = function(
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
) {
  expand.grid(
    sn          = sn_levels,
    acc_level   = accuracy,
    sig_level   = sigma_levels,
    omg_level   = omega_levels,
    mu_prev     = mu_prev_vals,
    tau         = tau,
    rho         = rho,
    sigma_bias  = sigma_bias,
    kappa_prior = kappa_prior,
    delta_se    = delta_se_vals,
    delta_sp    = delta_sp_vals
  )
}

# Scenario grid for Simulation 5
# prior_se / prior_sp: the single Se/Sp value used as the prior centre for
# ALL measures (regardless of their true per-measure accuracy)
scenario_grid_sim5 = function(
    sn_levels       = 20,
    accuracy        = c(0.80, 0.90),
    sigma_levels    = 0.2,
    omega_levels    = 0.2,
    mu_prev_vals    = c(0.015, 0.06, 0.12),
    tau             = 1.0,
    rho             = -0.5,
    sigma_bias      = 0.3,
    kappa_prior     = 200
) {
  # accuracy defines the common prior centres:
  #   prior_se = accuracy,  prior_sp = accuracy - 0.05
  expand.grid(
    sn          = sn_levels,
    acc_level   = accuracy,
    sig_level   = sigma_levels,
    omg_level   = omega_levels,
    mu_prev     = mu_prev_vals,
    tau         = tau,
    rho         = rho,
    sigma_bias  = sigma_bias,
    kappa_prior = kappa_prior
  )
}


# ---- Extended Dataset Generator -------------------------------------------

# Extends simulate_dataset to carry delta_se, delta_sp, prior_se, prior_sp
# through the data frame (these do NOT affect data generation, only fitting)
simulate_dataset_extended = function(
    sn = 20,
    meanlog_n = 5.50, sdlog_n = 0.55, add_const_n = 30,
    di = 10,
    mean_se = 0.80, mean_sp = 0.75,
    sd_se_logit = 0.1, sd_sp_logit = 0.1, rho = -0.5,
    omega_se = 0.2, omega_sp = 0.2,
    mu_prevalence = 0.015,
    tau_logit = 0.5,
    sigma_bias = 0, kappa_prior = 200,
    floor_prob = 0.5, ceiling_prob = 0.999,
    delta_se = 0, delta_sp = 0,
    prior_se = NA, prior_sp = NA
) {
  dat = simulate_dataset(
    sn = sn,
    meanlog_n = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
    di = di,
    mean_se = mean_se, mean_sp = mean_sp,
    sd_se_logit = sd_se_logit, sd_sp_logit = sd_sp_logit, rho = rho,
    omega_se = omega_se, omega_sp = omega_sp,
    mu_prevalence = mu_prevalence,
    tau_logit = tau_logit,
    sigma_bias = sigma_bias, kappa_prior = kappa_prior,
    floor_prob = floor_prob, ceiling_prob = ceiling_prob
  )
  dat$delta_se = delta_se
  dat$delta_sp = delta_sp
  dat$prior_se = prior_se
  dat$prior_sp = prior_sp
  dat
}


# ---- Simulation Generators (Sim 4 & 5) ------------------------------------

simulate_all_scenarios_sim4 = function(R = 200,
                                       base_seed = 20250903,
                                       meanlog_n = 5.50, sdlog_n = 0.55,
                                       add_const_n = 30,
                                       di = 10,
                                       floor_prob = 0.5,
                                       ceiling_prob = 0.999,
                                       grid = scenario_grid_sim4()) {
  n_scen = nrow(grid)
  out = vector("list", length = n_scen)

  for (s in seq_len(n_scen)) {
    message(sprintf("Scenario %d / %d", s, n_scen))

    sn       = grid$sn[s]
    acctag   = grid$acc_level[s]
    sigtag   = grid$sig_level[s]
    omgtag   = grid$omg_level[s]
    mu_p     = grid$mu_prev[s]
    tau_logit = grid$tau[s]
    c_rho    = grid$rho[s]
    c_sigma  = grid$sigma_bias[s]
    c_kappa  = grid$kappa_prior[s]
    c_dse    = grid$delta_se[s]
    c_dsp    = grid$delta_sp[s]

    mu_se = acctag
    mu_sp = acctag - 0.05
    sd_se = sigtag
    sd_sp = sigtag
    om_se = omgtag
    om_sp = omgtag

    datasets = vector("list", length = R)
    for (r in seq_len(R)) {
      set.seed(base_seed + (s - 1) * 1e6 + r)
      datasets[[r]] = simulate_dataset_extended(
        sn = sn,
        meanlog_n = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
        di = di,
        mean_se = mu_se, mean_sp = mu_sp,
        sd_se_logit = sd_se, sd_sp_logit = sd_sp, rho = c_rho,
        omega_se = om_se, omega_sp = om_sp,
        mu_prevalence = mu_p,
        tau_logit = tau_logit,
        sigma_bias = c_sigma, kappa_prior = c_kappa,
        floor_prob = floor_prob, ceiling_prob = ceiling_prob,
        delta_se = c_dse, delta_sp = c_dsp
      )
    }

    out[[s]] = list(
      scenario_id = s,
      factors = list(
        sn = sn, accuracy = acctag, sigma_level = sigtag,
        omega_level = omgtag, mu_prev = mu_p, rho = c_rho,
        delta_se = c_dse, delta_sp = c_dsp,
        sigma_bias = c_sigma, kappa_prior = c_kappa
      ),
      replicates = datasets
    )
  }
  out
}

simulate_all_scenarios_sim5 = function(R = 200,
                                       base_seed = 20250903,
                                       meanlog_n = 5.50, sdlog_n = 0.55,
                                       add_const_n = 30,
                                       di = 10,
                                       floor_prob = 0.5,
                                       ceiling_prob = 0.999,
                                       grid = scenario_grid_sim5()) {
  n_scen = nrow(grid)
  out = vector("list", length = n_scen)

  for (s in seq_len(n_scen)) {
    message(sprintf("Scenario %d / %d", s, n_scen))

    sn        = grid$sn[s]
    acctag    = grid$acc_level[s]
    sigtag    = grid$sig_level[s]
    omgtag    = grid$omg_level[s]
    mu_p      = grid$mu_prev[s]
    tau_logit = grid$tau[s]
    c_rho     = grid$rho[s]
    c_sigma   = grid$sigma_bias[s]
    c_kappa   = grid$kappa_prior[s]

    # The common prior centres for ALL measures
    prior_se = acctag
    prior_sp = acctag - 0.05

    # Data are generated with the SAME accuracy as the prior centre
    # (the bias arises because per-measure accuracy varies around this)
    mu_se = acctag
    mu_sp = acctag - 0.05
    sd_se = sigtag
    sd_sp = sigtag
    om_se = omgtag
    om_sp = omgtag

    datasets = vector("list", length = R)
    for (r in seq_len(R)) {
      set.seed(base_seed + (s - 1) * 1e6 + r)
      datasets[[r]] = simulate_dataset_extended(
        sn = sn,
        meanlog_n = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
        di = di,
        mean_se = mu_se, mean_sp = mu_sp,
        sd_se_logit = sd_se, sd_sp_logit = sd_sp, rho = c_rho,
        omega_se = om_se, omega_sp = om_sp,
        mu_prevalence = mu_p,
        tau_logit = tau_logit,
        sigma_bias = c_sigma, kappa_prior = c_kappa,
        floor_prob = floor_prob, ceiling_prob = ceiling_prob,
        prior_se = prior_se, prior_sp = prior_sp
      )
    }

    out[[s]] = list(
      scenario_id = s,
      factors = list(
        sn = sn, accuracy = acctag, sigma_level = sigtag,
        omega_level = omgtag, mu_prev = mu_p, rho = c_rho,
        prior_se = prior_se, prior_sp = prior_sp,
        sigma_bias = c_sigma, kappa_prior = c_kappa
      ),
      replicates = datasets
    )
  }
  out
}


# ---- Prior Builders ---------------------------------------------------

# Build M6 priors with SYSTEMATIC bias (Simulations 4a and b)
# Shifts each measure's prior centre by a fixed delta_se / delta_sp
# (probability scale) BEFORE applying the usual random perturbation
# (sigma_bias).
build_priors_m6_systematic = function(dat_fac,
                                      priors_common,
                                      n_eff        = 200,
                                      sigma_bias   = 0.3,
                                      rho_bias     = -0.5,
                                      delta_se     = 0,
                                      delta_sp     = 0,
                                      floor_prob   = 0.5,
                                      ceiling_prob = 0.999) {

  pri = c(priors_common$pi, priors_common$se_sp_sd)

  centers = dat_fac %>%
    dplyr::select(measure_id, se_base, sp_base) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  n_meas = nrow(centers)

  # Apply systematic bias on probability scale
  se_sys = centers$se_base + delta_se
  sp_sys = centers$sp_base + delta_sp

  # Clamp after systematic shift
  se_sys = clamp_probability(se_sys, floor_prob, ceiling_prob)
  sp_sys = clamp_probability(sp_sys, floor_prob, ceiling_prob)

  # Apply random bias on logit scale (as in build_priors_m6)
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_meas, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_meas), rep(0, n_meas))
  }

  eta_se = stats::qlogis(se_sys) + E[, 1]
  eta_sp = stats::qlogis(sp_sys) + E[, 2]
  se_b   = stats::plogis(eta_se)
  sp_b   = stats::plogis(eta_sp)

  # Final clamping
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

  for (i in seq_len(n_meas)) {
    m = centers$measure_id[i]
    pri = c(
      pri,
      prior_logit_for_measure(se_b[i], m, nlpar = "Se", n = n_eff),
      prior_logit_for_measure(sp_b[i], m, nlpar = "Sp", n = n_eff)
    )
  }

  pri
}

# Build M6 priors with a COMMON prior centre for all measures (Simulation 5a and b)
# Instead of using per-measure se_base/sp_base, centres ALL measure priors
# on a single prior_se / prior_sp.  Random perturbation (sigma_bias) is
# still applied per-measure in some cases.
build_priors_m6_common_center = function(dat_fac,
                                         priors_common,
                                         n_eff        = 200,
                                         sigma_bias   = 0.3,
                                         rho_bias     = -0.5,
                                         prior_se     = 0.85,
                                         prior_sp     = 0.80,
                                         floor_prob   = 0.5,
                                         ceiling_prob = 0.999) {

  pri = c(priors_common$pi, priors_common$se_sp_sd)

  centers = dat_fac %>%
    dplyr::select(measure_id) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  n_meas = nrow(centers)

  # Common centre for all measures (clamped)
  se_common = clamp_probability(prior_se, floor_prob, ceiling_prob)
  sp_common = clamp_probability(prior_sp, floor_prob, ceiling_prob)

  # Random bias per measure on logit scale
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_meas, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_meas), rep(0, n_meas))
  }

  eta_se = stats::qlogis(se_common) + E[, 1]
  eta_sp = stats::qlogis(sp_common) + E[, 2]
  se_b   = stats::plogis(eta_se)
  sp_b   = stats::plogis(eta_sp)

  # Final clamping
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

  for (i in seq_len(n_meas)) {
    m = centers$measure_id[i]
    pri = c(
      pri,
      prior_logit_for_measure(se_b[i], m, nlpar = "Se", n = n_eff),
      prior_logit_for_measure(sp_b[i], m, nlpar = "Sp", n = n_eff)
    )
  }

  pri
}


# ---- Fitting Functions (M1 & M6 only) ------------------------------------

# Fit M1 and M6 for one replicate - Simulation 4
# M6 priors incorporate systematic bias (delta_se, delta_sp) plus random
# perturbation (sigma_bias). I decided to break this into 4a and 4b in the end.
#  I didn't bother implementing bounded variants for these simulations.
fit_m1_m6_sim4 = function(dat_raw,
                          iter = 6000,
                          adapt_delta = 0.999,
                          max_treedepth = 20,
                          stepsize = NULL) {

  true_pr   = dat_raw$mu_prevalence[1]
  kappa_pr  = dat_raw$kappa_prior[1]
  c_sigma   = dat_raw$sigma_bias[1]
  c_delta_se = dat_raw$delta_se[1]
  c_delta_sp = dat_raw$delta_sp[1]

  common_priors = make_common_priors(qlogis(true_pr))

  if (is.null(stepsize)) {
    ctrl = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  } else {
    ctrl = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                step_size = stepsize)
  }

  # M1: Naive model
  dat_m1 = dat_raw %>% dplyr::select(y, n, study)

  fit1 = brm(
    y | trials(n) ~ 1 + (1 | study),
    data    = dat_m1,
    family  = binomial(),
    backend = 'cmdstanr',
    prior   = common_priors$m123,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # M6: Misclassification model with systematic + random bias
  dat_m6 = dat_raw %>% dplyr::select(y, n, study, measure_id, se_base, sp_base)

  pri6 = build_priors_m6_systematic(
    dat_m6, common_priors,
    n_eff      = kappa_pr,
    sigma_bias = c_sigma,
    delta_se   = c_delta_se,
    delta_sp   = c_delta_sp
  )

  fit6 = brm(
    formula = bf_m6,
    data    = dat_m6 %>% dplyr::mutate(measure_id = as.character(measure_id)),
    family  = binomial(link = 'identity'),
    backend = 'cmdstanr',
    prior   = pri6,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # Collect metrics (I'm still calling this matrix k, which was from when I was
  # prototyping my first sims, sadly it's stuck at this stage!)
  k = dplyr::bind_rows(
    simulation_metrics(fit1, dat_raw) %>% dplyr::mutate(model = "M1"),
    simulation_metrics(fit6, dat_raw) %>% dplyr::mutate(model = "M6")
  ) %>%
    dplyr::bind_cols(
      dat_raw %>%
        dplyr::select(sn, di, mean_se, mean_sp, sd_se_logit, sd_sp_logit,
                      rho, omega_se, omega_sp, sigma_bias, kappa_prior,
                      mu_prevalence, tau_logit, delta_se, delta_sp) %>%
        dplyr::distinct()
    )

  return(k)
}


# Fit M1 and M6 for one replicate - Simulation 5
# M6 priors use a COMMON centre (prior_se / prior_sp) for all measures,
# with or without random perturbation (sigma_bias). I decided
# to break this into 5a and 5b in the end. I didn't bother
# implementing bounded variants for these simulations.
fit_m1_m6_sim5 = function(dat_raw,
                          iter = 6000,
                          adapt_delta = 0.999,
                          max_treedepth = 20,
                          stepsize = NULL) {

  true_pr   = dat_raw$mu_prevalence[1]
  kappa_pr  = dat_raw$kappa_prior[1]
  c_sigma   = dat_raw$sigma_bias[1]
  c_prior_se = dat_raw$prior_se[1]
  c_prior_sp = dat_raw$prior_sp[1]

  common_priors = make_common_priors(qlogis(true_pr))

  if (is.null(stepsize)) {
    ctrl = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  } else {
    ctrl = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                step_size = stepsize)
  }

  # M1: Naive model
  dat_m1 = dat_raw %>% dplyr::select(y, n, study)

  fit1 = brm(
    y | trials(n) ~ 1 + (1 | study),
    data    = dat_m1,
    family  = binomial(),
    backend = 'cmdstanr',
    prior   = common_priors$m123,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # M6: Misclassification model with common-centre priors
  dat_m6 = dat_raw %>% dplyr::select(y, n, study, measure_id, se_base, sp_base)

  pri6 = build_priors_m6_common_center(
    dat_m6, common_priors,
    n_eff      = kappa_pr,
    sigma_bias = c_sigma,
    prior_se   = c_prior_se,
    prior_sp   = c_prior_sp
  )

  fit6 = brm(
    formula = bf_m6,
    data    = dat_m6 %>% dplyr::mutate(measure_id = as.character(measure_id)),
    family  = binomial(link = 'identity'),
    backend = 'cmdstanr',
    prior   = pri6,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # Collect metrics (again, apologies for bad naming convention)
  k = dplyr::bind_rows(
    simulation_metrics(fit1, dat_raw) %>% dplyr::mutate(model = "M1"),
    simulation_metrics(fit6, dat_raw) %>% dplyr::mutate(model = "M6")
  ) %>%
    dplyr::bind_cols(
      dat_raw %>%
        dplyr::select(sn, di, mean_se, mean_sp, sd_se_logit, sd_sp_logit,
                      rho, omega_se, omega_sp, sigma_bias, kappa_prior,
                      mu_prevalence, tau_logit, prior_se, prior_sp) %>%
        dplyr::distinct()
    )

  return(k)
}


# ---- Run / Repair Functions (M1 & M6 only) --------------------------------

# Used to fit the simulations that use only M1 or M6 (not the
# other models); this should be Simulations 4 - 6. Here I figured out
# I could pass the function for each, so it is accepted as an argument
# instead of making two functions for handle the different sims.
run_all_rds_m1m6 = function(sims, dir, fit_fn,
                            iter = 6000, adapt_delta = 0.999,
                            skip_num = 0) {
  for (i in (seq_len(length(sims) - skip_num) + skip_num)) {
    message(sprintf("Starting scenario i = %d", i))
    c_sim = sims[[i]]

    for (j in seq_len(length(c_sim$replicates))) {
      fn = sprintf('%sSim%d_%d.rds', dir, i, j)

      if (file.exists(fn)) {
        message(sprintf("  Loading r = %d", j)) # Note it says loading because my old code did load, legacy now, it skips
        next
      }

      start_time = Sys.time()
      c_rep = c_sim$replicates[[j]]

      summ = fit_fn(c_rep, iter = iter, adapt_delta = adapt_delta)

      temp = summ %>%
        dplyr::mutate(sim_id = i, rep_id = j) %>%
        dplyr::select(sim_id, rep_id, model, everything()) %>%
        dplyr::mutate(data = list(c_rep))

      saveRDS(temp, fn)
      end_time = Sys.time()
      message(sprintf("     Completed r = %d in %.02f s", j,
                      as.numeric(end_time - start_time, units = "secs")))
    }
  }
}

# Repair cycle for M1 + M6 simulations, see the original repair function
# for further details.
repair_cycle_m1m6 = function(sims, fit_fn,
                             iter = 12000, adapt_delta = 0.999,
                             step_size = 0.002,
                             dir = 'simulation_data/',
                             cor_dir = 'corrected_data/',
                             cycles = 3) {
  for (cyc in seq_len(cycles)) {
    message(sprintf("=== Repair cycle %d ===", cyc))
    results_dat = read_sim_rds(dir = dir)

    bad_trials = results_dat %>%
      dplyr::filter(ndiv > 0 | maxrhat > 1.01) %>%
      dplyr::select(sim_id, rep_id) %>%
      dplyr::distinct()

    # Clear old corrected files
    for (j in seq_len(nrow(bad_trials))) {
      bad_fn = sprintf('%sSim%d_%d.rds', cor_dir,
                       bad_trials$sim_id[j], bad_trials$rep_id[j])
      if (file.exists(bad_fn)) file.remove(bad_fn)
    }

    if (nrow(bad_trials) > 0) {
      fix_all_rds_m1m6(bad_trials, sims, fit_fn,
                       iter = iter, adapt_delta = adapt_delta,
                       step_size = step_size,
                       dir = dir, cor_dir = cor_dir)
    } else {
      return(0)
    }
  }

  # Final count of remaining bad trials
  results_dat = read_sim_rds(dir = dir)
  bad_trials = results_dat %>%
    dplyr::filter(ndiv > 0 | maxrhat > 1.01) %>%
    dplyr::select(sim_id, rep_id) %>%
    dplyr::distinct()
  return(nrow(bad_trials))
}

fix_all_rds_m1m6 = function(bad_trials, sims, fit_fn,
                            iter = 12000, adapt_delta = 0.999,
                            step_size = 0.002,
                            dir = 'simulation_data/',
                            cor_dir = 'corrected_data/') {
  # Keep the previous result before replacing it with a repaired fit.
  if (!dir.exists(cor_dir) && !dir.create(cor_dir, recursive = TRUE)) {
    stop("Could not create backup directory: ", cor_dir)
  }

  nbad = nrow(bad_trials)
  for (i in seq_len(nbad)) {
    s = bad_trials$sim_id[i]
    r = bad_trials$rep_id[i]
    message(sprintf("Fixing s=%d, r=%d  (%d of %d)", s, r, i, nbad))

    fn     = sprintf('%sSim%d_%d.rds', dir, s, r)
    bad_fn = sprintf('%sSim%d_%d.rds', cor_dir, s, r)

    if (file.exists(bad_fn)) {
      message("    File exists in corrected dir. Skipping.")
      next
    }

    if (!file.rename(fn, bad_fn)) {
      stop("Could not back up ", fn, "; the original result has not been replaced.")
    }

    dat_raw = sims[[s]]$replicates[[r]]

    start_time = Sys.time()
    summ = fit_fn(dat_raw, iter = iter, adapt_delta = adapt_delta,
                  stepsize = step_size)

    temp = summ %>%
      dplyr::mutate(sim_id = s, rep_id = r) %>%
      dplyr::select(sim_id, rep_id, model, everything()) %>%
      dplyr::mutate(data = list(dat_raw))

    saveRDS(temp, fn)
    end_time = Sys.time()
    message(sprintf("     Completed s=%d, r=%d in %.2f s", s, r,
                    as.numeric(end_time - start_time, units = "secs")))
  }
}