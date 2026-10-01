# ============================================================================
# Simulation 6: Mixed Gold Standard and Screening Studies
# ============================================================================
#
# Evaluates whether the inclusion of gold standard studies (Se = Sp ≈ 1)
# can rescue prevalence estimation under systematic Sp misspecification.
#
# Design: 3 (proportion gold: .10, .25, .50) ×
#         3 (true prevalence: .015, .06, .12) ×
#         3 (δSp: −0.10, 0, +0.10) = 27 scenarios
#
# δSe is fixed at 0 (Sim 4 showed Se misspecification is secondary).
# Only M1 and M6 are fit.
# ============================================================================


# ---- Scenario Grid --------------------------------------------------------

scenario_grid_sim6 = function(
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
    delta_sp    = delta_sp_vals,
    prop_gold   = prop_gold_vals
  )
}


# ---- Data Generator with Gold Standard Studies ----------------------------

# Generate a dataset containing a mix of gold standard and screening studies.
# Gold standard studies have Se = Sp = ceiling_prob (effectively perfect
# classification). They are assigned a special measure_id = 0.
# Screening studies are generated as usual via simulate_dataset().
simulate_dataset_gold = function(
    sn = 20,
    meanlog_n = 5.50, sdlog_n = 0.55, add_const_n = 30,
    di = 10,
    mean_se = 0.85, mean_sp = 0.80,
    sd_se_logit = 0.2, sd_sp_logit = 0.2, rho = -0.5,
    omega_se = 0.2, omega_sp = 0.2,
    mu_prevalence = 0.015,
    tau_logit = 1.0,
    sigma_bias = 0.3, kappa_prior = 200,
    floor_prob = 0.5, ceiling_prob = 0.999,
    prop_gold = 0.10,
    delta_sp = 0
) {
  # Determine number of gold vs screening studies
  n_gold     = max(1, round(sn * prop_gold))   # at least 1 gold study
  n_screen   = sn - n_gold

  # Gold standard studies
  n_gold_i   = sample_sizes_logn(n_gold, meanlog_n, sdlog_n, add_const_n)
  theta_gold = simulate_true_prevalence(n_gold, mu_prevalence, tau_logit)

  # Gold standard: no misclassification (Se = Sp = near-perfect)
  Se_gold = rep(ceiling_prob, n_gold)
  Sp_gold = rep(ceiling_prob, n_gold)

  fm_gold = forward_misclassification(n_gold_i, theta_gold, Se_gold, Sp_gold)

  dat_gold = data.frame(
    sn            = sn,
    di            = di,
    mean_se       = mean_se,
    mean_sp       = mean_sp,
    sd_se_logit   = sd_se_logit,
    sd_sp_logit   = sd_sp_logit,
    rho           = rho,
    omega_se      = omega_se,
    omega_sp      = omega_sp,
    sigma_bias    = sigma_bias,
    kappa_prior   = kappa_prior,
    mu_prevalence = mu_prevalence,
    tau_logit     = tau_logit,
    study         = seq_len(n_gold),
    n             = n_gold_i,
    y             = fm_gold$y,
    p_obs         = fm_gold$p_obs,
    theta_true    = theta_gold,
    Se            = Se_gold,
    Sp            = Sp_gold,
    measure_id    = 0,            # special ID for gold standard
    se_base       = ceiling_prob,
    sp_base       = ceiling_prob,
    eta_se_base   = qlogis(ceiling_prob),
    eta_sp_base   = qlogis(ceiling_prob),
    is_gold       = TRUE,
    prop_gold     = prop_gold,
    delta_sp      = delta_sp,
    stringsAsFactors = FALSE
  )

  # Screening studies
  if (n_screen > 0) {
    # Generate screening data using existing function
    dat_screen = simulate_dataset(
      sn = n_screen,
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

    # Shift study IDs to follow gold studies
    dat_screen$study = dat_screen$study + n_gold
    dat_screen$is_gold  = FALSE
    dat_screen$prop_gold = prop_gold
    dat_screen$delta_sp  = delta_sp
    dat_screen$sn = sn # Fixes a multiple row issue owing to SN mismatch


    dat = rbind(dat_gold, dat_screen)
  } else {
    dat = dat_gold
  }

  # Re-index studies
  dat$study = seq_len(nrow(dat))
  dat
}

# ---- Simulation Generator -------------------------------------------------

simulate_all_scenarios_sim6 = function(R = 200,
                                       base_seed = 20250903,
                                       meanlog_n = 5.50, sdlog_n = 0.55,
                                       add_const_n = 30,
                                       di = 10,
                                       floor_prob = 0.5,
                                       ceiling_prob = 0.999,
                                       grid = simulation6) {
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
    c_dsp     = grid$delta_sp[s]
    c_pgold   = grid$prop_gold[s]

    mu_se = acctag
    mu_sp = acctag - 0.05

    datasets = vector("list", length = R)
    for (r in seq_len(R)) {
      set.seed(base_seed + (s - 1) * 1e6 + r)
      datasets[[r]] = simulate_dataset_gold(
        sn = sn,
        meanlog_n = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
        di = di,
        mean_se = mu_se, mean_sp = mu_sp,
        sd_se_logit = sigtag, sd_sp_logit = sigtag, rho = c_rho,
        omega_se = omgtag, omega_sp = omgtag,
        mu_prevalence = mu_p,
        tau_logit = tau_logit,
        sigma_bias = c_sigma, kappa_prior = c_kappa,
        floor_prob = floor_prob, ceiling_prob = ceiling_prob,
        prop_gold = c_pgold,
        delta_sp = c_dsp
      )
    }

    out[[s]] = list(
      scenario_id = s,
      factors = list(
        sn = sn, accuracy = acctag, sigma_level = sigtag,
        omega_level = omgtag, mu_prev = mu_p, rho = c_rho,
        delta_sp = c_dsp, prop_gold = c_pgold,
        sigma_bias = c_sigma, kappa_prior = c_kappa
      ),
      replicates = datasets
    )
  }
  out
}


# ---- Prior Builder for Sim 6 ----------------------------------------------

# Build M6 priors for mixed gold/screening data
# Screening measures: priors shifted by delta_sp (prob scale) then
#   perturbed by sigma_bias (logit scale), as in build_priors_m6_systematic.
# Gold measures (measure_id == 0): priors fixed at constant(10) on logit
#   scale, matching the case study approach.

build_priors_m6_gold = function(dat_fac,
                                priors_common,
                                n_eff        = 200,
                                sigma_bias   = 0.3,
                                rho_bias     = -0.5,
                                delta_sp     = 0,
                                floor_prob   = 0.5,
                                ceiling_prob = 0.999) {

  pri = c(priors_common$pi, priors_common$se_sp_sd)

  centers = dat_fac %>%
    dplyr::select(measure_id, se_base, sp_base, is_gold) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  # Separate gold and screening measures
  gold_centers   = centers %>% dplyr::filter(is_gold)
  screen_centers = centers %>% dplyr::filter(!is_gold)

  # Gold standard priors: fixed at near-perfect
  for (i in seq_len(nrow(gold_centers))) {
    m = gold_centers$measure_id[i]
    pri = c(
      pri,
      set_prior("constant(10)", nlpar = "Se", class = "b",
                coef = sprintf("measure_id%d", m)),
      set_prior("constant(10)", nlpar = "Sp", class = "b",
                coef = sprintf("measure_id%d", m))
    )
  }

  # Screening measure priors: systematic + random bias
  n_screen = nrow(screen_centers)

  if (n_screen > 0) {
    # Systematic Sp shift (no Se shift: delta_se = 0)
    se_sys = screen_centers$se_base
    sp_sys = screen_centers$sp_base + delta_sp

    se_sys = clamp_probability(se_sys, floor_prob, ceiling_prob)
    sp_sys = clamp_probability(sp_sys, floor_prob, ceiling_prob)

    # Random perturbation
    if (sigma_bias > 0) {
      Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
      E = MASS::mvrnorm(n = n_screen, mu = c(0, 0), Sigma = Sigma)
    } else {
      E = cbind(rep(0, n_screen), rep(0, n_screen))
    }

    eta_se = stats::qlogis(se_sys) + E[, 1]
    eta_sp = stats::qlogis(sp_sys) + E[, 2]
    se_b   = stats::plogis(eta_se)
    sp_b   = stats::plogis(eta_sp)

    se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
    sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

    for (i in seq_len(n_screen)) {
      m = screen_centers$measure_id[i]
      pri = c(
        pri,
        prior_logit_for_measure(se_b[i], m, nlpar = "Se", n = n_eff),
        prior_logit_for_measure(sp_b[i], m, nlpar = "Sp", n = n_eff)
      )
    }
  }

  pri
}


# ---- Fitting Function (M1 & M6) -------------------------------------------

fit_m1_m6_sim6 = function(dat_raw,
                          iter = 6000,
                          adapt_delta = 0.999,
                          max_treedepth = 20,
                          stepsize = NULL) {

  true_pr   = dat_raw$mu_prevalence[1]
  kappa_pr  = dat_raw$kappa_prior[1]
  c_sigma   = dat_raw$sigma_bias[1]
  c_dsp     = dat_raw$delta_sp[1]

  common_priors = make_common_priors(qlogis(true_pr))

  if (is.null(stepsize)) {
    ctrl = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  } else {
    ctrl = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                step_size = stepsize)
  }

  # M1: Naive model (ignores gold/screening distinction)
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

  # M6: Misclassification model with gold + screening
  dat_m6 = dat_raw %>%
    dplyr::select(y, n, study, measure_id, se_base, sp_base, is_gold)

  pri6 = build_priors_m6_gold(
    dat_m6, common_priors,
    n_eff      = kappa_pr,
    sigma_bias = c_sigma,
    delta_sp   = c_dsp
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

  # Collect metrics
  dat_raw$sn = dat_raw$sn[1] # Fixes a multiple row issue owing to SN mismatch

  k = dplyr::bind_rows( # Again, sorry for k, I may fix it for the next one
    simulation_metrics(fit1, dat_raw)  %>% dplyr::mutate(model = "M1"),
    simulation_metrics(fit6, dat_raw)  %>% dplyr::mutate(model = "M6")
  ) %>%
    dplyr::bind_cols(
      dat_raw %>%
        dplyr::select(sn, di, mean_se, mean_sp, sd_se_logit, sd_sp_logit,
                      rho, omega_se, omega_sp, sigma_bias, kappa_prior,
                      mu_prevalence, tau_logit, delta_sp, prop_gold) %>%
        dplyr::distinct() %>% dplyr::slice(1)
    )

  return(k)
}
