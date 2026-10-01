# ============================================================================
# Simulation 8: Se/Sp Estimation and Recovery
# ============================================================================
#
# Evaluates two approaches for estimating Se/Sp from the data when
# gold-standard studies are available:
#   1. Discrete KWGA (kernel-weighted grid averaging): fast algorithm scoring
#      Se/Sp pairs by gold-screener agreement in an uncorrected comparison model
#   2. Joint model: brms model estimating a common Se/Sp for screeners
#      with gold studies anchoring prevalence
#
# For each approach, the estimated Se/Sp is used to set priors for M6,
# and prevalence recovery is evaluated.
#
# Design: 3 (Sp: .70, .80, .90) ×
#         3 (prevalence: .015, .06, .12) ×
#         3 (prop gold: .10, .25, .50) ×
#         2 (k: 20, 40) = 54 scenarios
#
# Se fixed at .85 throughout. sigma_bias = 0, delta_sp = 0.
# ============================================================================


# ---- Scenario Grid --------------------------------------------------------

scenario_grid_sim8 = function(
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
) {
  expand.grid(
    sn          = sn_levels,
    mean_se     = mean_se_vals,
    mean_sp     = mean_sp_vals,
    sig_level   = sigma_levels,
    omg_level   = omega_levels,
    mu_prev     = mu_prev_vals,
    tau         = tau,
    rho         = rho,
    sigma_bias  = sigma_bias,
    kappa_prior = kappa_prior,
    prop_gold   = prop_gold_vals,
    stringsAsFactors = FALSE
  )
}


# ---- Data Generator (reuses simulate_dataset_gold from Sim 6) -------------

# Wrapper for simulate_dataset_gold that accepts separate Se and Sp
# rather than the accuracy/accuracy-0.05 convention from Sim 6.

simulate_dataset_sim8 = function(
    sn = 20,
    meanlog_n = 5.50, sdlog_n = 0.55, add_const_n = 30,
    di = 10,
    mean_se = 0.85, mean_sp = 0.80,
    sd_se_logit = 0.2, sd_sp_logit = 0.2, rho = -0.5,
    omega_se = 0.2, omega_sp = 0.2,
    mu_prevalence = 0.015,
    tau_logit = 1.0,
    sigma_bias = 0, kappa_prior = 100,
    floor_prob = 0.5, ceiling_prob = 0.999,
    prop_gold = 0.10
) {
  simulate_dataset_gold(
    sn = sn,
    meanlog_n = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
    di = di,
    mean_se = mean_se, mean_sp = mean_sp,
    sd_se_logit = sd_se_logit, sd_sp_logit = sd_sp_logit, rho = rho,
    omega_se = omega_se, omega_sp = omega_sp,
    mu_prevalence = mu_prevalence,
    tau_logit = tau_logit,
    sigma_bias = sigma_bias, kappa_prior = kappa_prior,
    floor_prob = floor_prob, ceiling_prob = ceiling_prob,
    prop_gold = prop_gold,
    delta_sp = 0
  )
}


# ---- Simulation Generator -------------------------------------------------

simulate_all_scenarios_sim8 = function(R = 200,
                                        base_seed = 20250321,
                                        meanlog_n = 5.50, sdlog_n = 0.55,
                                        add_const_n = 30,
                                        di = 10,
                                        floor_prob = 0.5,
                                        ceiling_prob = 0.999,
                                        grid = NULL) {
  n_scen = nrow(grid)
  out = vector("list", length = n_scen)

  for (s in seq_len(n_scen)) {
    message(sprintf("Scenario %d / %d", s, n_scen))

    datasets = vector("list", length = R)
    for (r in seq_len(R)) {
      set.seed(base_seed + (s - 1) * 1e6 + r)
      datasets[[r]] = simulate_dataset_sim8(
        sn            = grid$sn[s],
        meanlog_n     = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
        di            = di,
        mean_se       = grid$mean_se[s],
        mean_sp       = grid$mean_sp[s],
        sd_se_logit   = grid$sig_level[s],
        sd_sp_logit   = grid$sig_level[s],
        rho           = grid$rho[s],
        omega_se      = grid$omg_level[s],
        omega_sp      = grid$omg_level[s],
        mu_prevalence = grid$mu_prev[s],
        tau_logit     = grid$tau[s],
        sigma_bias    = grid$sigma_bias[s],
        kappa_prior   = grid$kappa_prior[s],
        floor_prob    = floor_prob,
        ceiling_prob  = ceiling_prob,
        prop_gold     = grid$prop_gold[s]
      )
    }

    out[[s]] = list(
      scenario_id = s,
      factors = list(
        sn = grid$sn[s],
        mean_se = grid$mean_se[s],
        mean_sp = grid$mean_sp[s],
        sigma_level = grid$sig_level[s],
        omega_level = grid$omg_level[s],
        mu_prev = grid$mu_prev[s],
        tau = grid$tau[s],
        rho = grid$rho[s],
        prop_gold = grid$prop_gold[s],
        sigma_bias = grid$sigma_bias[s],
        kappa_prior = grid$kappa_prior[s]
      ),
      replicates = datasets
    )
  }
  out
}


# ---- Fitting Functions: Sim 8 --------------------------------------------

# For each replicate, fit_sim8() fits the comparison, joint, M6d and M6j
# models, computes the KWGA estimates, and returns one wide row:
#   1. Uncorrected comparison model (gold vs screener)       fit_comparison_sim8()
#   2. Discrete KWGA estimates of Se/Sp and prevalence        estimate_sesp_kwga_sim8()
#   3. Joint brms model estimating Se/Sp                      fit_joint_sim8()
#   4. M6 with priors centred on the KWGA estimates (M6d)     fit_m6_from_estimates_sim8()
#   5. M6 with priors centred on the joint estimates (M6j)    fit_m6_from_estimates_sim8()
# The summarise_*_sim8() helpers turn each fit into its block of result
# columns and sim8_result_columns fixes the column order of the stored files.
#
# seed: optional integer. When supplied, each model receives its own derived
# seed (seed + 1, ..., seed + 5) so that a replicate can be refitted
# reproducibly; NA (the default) leaves the samplers unseeded. The runners
# derive it per replicate from seed_base (see sim8_replicate_seed()), and each
# repair cycle adds cycle x 1e8 so that a failed repair is retried with a
# different, still reproducible, draw.

sim8_control = function(adapt_delta, max_treedepth, stepsize) {
  if (is.null(stepsize)) {
    list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  } else {
    list(adapt_delta = adapt_delta, max_treedepth = max_treedepth,
         step_size = stepsize)
  }
}

sim8_seed = function(seed, k) {
  if (is.na(seed)) NA else as.integer(seed + k)
}

# Uncorrected comparison model (gold vs screener); returns the brmsfit
fit_comparison_sim8 = function(dat_raw,
                               iter = 6000,
                               adapt_delta = 0.999,
                               max_treedepth = 20,
                               stepsize = NULL,
                               seed = NA) {

  true_pr = dat_raw$mu_prevalence[1]

  dat_comp = dat_raw %>%
    dplyr::select(y, n, study, is_gold) %>%
    dplyr::mutate(es_id = factor(study))

  true_app = sum(dat_comp$y[!dat_comp$is_gold]) / sum(dat_comp$n[!dat_comp$is_gold])      # screening-subgroup apparent rate based on observed
  # true_app = true_se * true_pr + (1 - true_sp) * (1 - true_pr)   #apparent screening rate; could be an alternate for the uncorrected prevalence but assumes true values known

  brm(
    y | trials(n) ~ 0 + is_gold + (1 | es_id),
    data    = dat_comp,
    family  = binomial(),
    backend = "cmdstanr",
    prior   = c(
      set_prior(sprintf("normal(%f, 1.5)", safe_qlogis(true_app)),
                class = "b", coef = "is_goldFALSE"),
      set_prior(sprintf("normal(%f, 1.5)", qlogis(true_pr)),
                class = "b", coef = "is_goldTRUE"),
      prior(normal(0, 1), class = "sd")
    ),
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = sim8_control(adapt_delta, max_treedepth, stepsize),
    seed    = seed
  )
}

# Discrete KWGA: Se/Sp grid weights, point estimates and the KWGA prevalence;
# returns one row of disc_* and kwga_* columns
estimate_sesp_kwga_sim8 = function(fit_comp, dat_raw,
                                   se_grid = seq(0.60, 0.95, by = 0.025),
                                   sp_grid = seq(0.60, 0.95, by = 0.025),
                                   n_draws = 4000,
                                   seed = NA) {

  true_se = dat_raw$mean_se[1]
  true_sp = dat_raw$mean_sp[1]
  true_pr = dat_raw$mu_prevalence[1]

  weights = compute_sesp_weights(
    fit_comp = fit_comp,
    se_vals  = se_grid,
    sp_vals  = sp_grid
  )

  # Extract discrete estimates
  discrete_joint_map <- weights$grid %>% dplyr::slice(1)

  discrete_marginal_sp <- weights$grid %>%
    dplyr::group_by(sp) %>%
    dplyr::summarise(w = sum(weight), .groups = "drop") %>%
    dplyr::slice_max(w, n = 1)

  discrete_marginal_se <- weights$grid %>%
    dplyr::group_by(se) %>%
    dplyr::summarise(w = sum(weight), .groups = "drop") %>%
    dplyr::slice_max(w, n = 1)

  discrete_weighted_se <- sum(weights$grid$se * weights$grid$weight)
  discrete_weighted_sp <- sum(weights$grid$sp * weights$grid$weight)

  # KWGA prevalence - pass the grid, not the list:
  if (!is.na(seed)) set.seed(seed)
  kwga_result <- kwga_prevalence(weights$grid, fit_comp = fit_comp,
                                 mode = "analytic", n_draws = n_draws)

  tibble::tibble(
    # --- Discrete KWGA estimates ---
    disc_joint_map_se  = discrete_joint_map$se,
    disc_joint_map_sp  = discrete_joint_map$sp,
    disc_marginal_se   = discrete_marginal_se$se[1],
    disc_marginal_sp   = discrete_marginal_sp$sp[1],
    disc_weighted_se   = discrete_weighted_se,
    disc_weighted_sp   = discrete_weighted_sp,
    disc_top_weight    = discrete_joint_map$weight,

    # Discrete Se/Sp bias
    disc_se_bias_wt    = discrete_weighted_se - true_se,
    disc_sp_bias_wt    = discrete_weighted_sp - true_sp,
    disc_sp_bias_marg  = discrete_marginal_sp$sp[1] - true_sp,

    # KWGA prevalence
    kwga_prev_mean      = kwga_result$summary$mean,
    kwga_prev_lb        = kwga_result$summary$ci_lb,
    kwga_prev_ub        = kwga_result$summary$ci_ub,
    kwga_prev_bias      = kwga_result$summary$mean - true_pr,
    kwga_prev_coverage  = kwga_result$summary$ci_lb < true_pr &
      kwga_result$summary$ci_ub > true_pr,
    kwga_n_eff_models   = kwga_result$summary$n_effective_models
  )
}

# Sampler diagnostics of the comparison model
summarise_comparison_sim8 = function(fit_comp) {
  tibble::tibble(
    comp_maxrhat = max(rhat(fit_comp), na.rm = TRUE),
    comp_ndiv    = n_divergent(fit_comp)
  )
}

# Joint brms model estimating a common screening Se/Sp; returns the brmsfit
fit_joint_sim8 = function(dat_raw,
                          iter = 6000,
                          adapt_delta = 0.999,
                          max_treedepth = 20,
                          stepsize = NULL,
                          seed = NA) {

  true_pr = dat_raw$mu_prevalence[1]

  dat_joint = dat_raw %>%
    dplyr::select(y, n, study, is_gold) %>%
    dplyr::mutate(
      es_id = factor(study),
      acc_group = factor(ifelse(is_gold, "gold", "screen"))
    )

  # Map Se/Sp to bounded inner logit scale
  se_inner_prior = qlogis(map_to_inner_t(0.80))
  sp_inner_prior = qlogis(map_to_inner_t(0.80))

  brm(
    formula = brms::bf(
      y | trials(n) ~ inv_logit(pi) * (0.5 + 0.5 * inv_logit(Se)) +
        (1 - inv_logit(pi)) * (1 - (0.5 + 0.5 * inv_logit(Sp))),
      pi ~ 1 + (1 | s | es_id),
      Se ~ 0 + acc_group + (1 | s | es_id),
      Sp ~ 0 + acc_group + (1 | s | es_id),
      nl = TRUE
    ),
    data    = dat_joint,
    family  = binomial(link = "identity"),
    backend = "cmdstanr",
    prior   = c(
      set_prior(sprintf("normal(%f, 1.5)", qlogis(true_pr)), nlpar = "pi"),
      prior(normal(0, 1), nlpar = "pi", class = "sd"),
      prior(normal(0, 0.5), nlpar = "Se", class = "sd"),
      prior(normal(0, 0.5), nlpar = "Sp", class = "sd"),
      prior_string("constant(10)", class = "b",
                   coef = "acc_groupgold", nlpar = "Se"),
      prior_string("constant(10)", class = "b",
                   coef = "acc_groupgold", nlpar = "Sp"),
      prior_string(sprintf("normal(%f, 1.0)", se_inner_prior),
                   class = "b", coef = "acc_groupscreen", nlpar = "Se"),
      prior_string(sprintf("normal(%f, 1.0)", sp_inner_prior),
                   class = "b", coef = "acc_groupscreen", nlpar = "Sp")
    ),
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = sim8_control(adapt_delta, max_treedepth, stepsize),
    seed    = seed
  )
}

# Joint-model Se/Sp and prevalence recovery (bounded scale) plus diagnostics
summarise_joint_sim8 = function(fit_joint, dat_raw) {

  true_se = dat_raw$mean_se[1]
  true_sp = dat_raw$mean_sp[1]
  true_pr = dat_raw$mu_prevalence[1]

  # Extract joint estimates (bounded scale)
  joint_draws = as_draws_df(fit_joint)
  joint_nms = names(joint_draws)
  se_col = grep("Se.*acc_groupscreen", joint_nms, value = TRUE)[1]
  sp_col = grep("Sp.*acc_groupscreen", joint_nms, value = TRUE)[1]

  joint_se_draws = 0.5 + 0.5 * plogis(joint_draws[[se_col]])
  joint_sp_draws = 0.5 + 0.5 * plogis(joint_draws[[sp_col]])
  joint_prev_draws = plogis(joint_draws$b_pi_Intercept)

  tibble::tibble(
    joint_se_mean      = mean(joint_se_draws),
    joint_sp_mean      = mean(joint_sp_draws),
    joint_se_sd        = sd(joint_se_draws),
    joint_sp_sd        = sd(joint_sp_draws),
    joint_se_bias      = mean(joint_se_draws) - true_se,
    joint_sp_bias      = mean(joint_sp_draws) - true_sp,
    joint_se_coverage  = quantile(joint_se_draws, 0.025) < true_se &
      quantile(joint_se_draws, 0.975) > true_se,
    joint_sp_coverage  = quantile(joint_sp_draws, 0.025) < true_sp &
      quantile(joint_sp_draws, 0.975) > true_sp,
    joint_prev_mean    = mean(joint_prev_draws),
    joint_prev_bias    = mean(joint_prev_draws) - true_pr,
    joint_prev_coverage = quantile(joint_prev_draws, 0.025) < true_pr &
      quantile(joint_prev_draws, 0.975) > true_pr,
    joint_maxrhat      = max(rhat(fit_joint), na.rm = TRUE),
    joint_ndiv         = n_divergent(fit_joint)
  )
}

# M6 (bounded) with screening Se/Sp priors centred on estimated values; used
# for M6d (KWGA estimates) and M6j (joint-model estimates). Returns the brmsfit
fit_m6_from_estimates_sim8 = function(dat_raw, se_est, sp_est,
                                      iter = 6000,
                                      adapt_delta = 0.999,
                                      max_treedepth = 20,
                                      stepsize = NULL,
                                      seed = NA) {

  true_pr  = dat_raw$mu_prevalence[1]
  kappa_pr = dat_raw$kappa_prior[1]
  common_priors = make_common_priors(qlogis(true_pr))

  dat_m6 = dat_raw %>%
    dplyr::select(y, n, study, measure_id, se_base, sp_base, is_gold) %>%
    dplyr::mutate(measure_id = as.character(measure_id))

  # Override the screening Se/Sp class-level priors with the estimates
  pri_override = build_m6_priors_from_estimates_bounded(
    dat_raw, common_priors,
    se_est = se_est,
    sp_est = sp_est,
    kappa = kappa_pr
  )

  brm(
    formula = bf_m6_bounded,
    data    = dat_m6,
    family  = binomial(link = "identity"),
    backend = "cmdstanr",
    prior   = pri_override,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = sim8_control(adapt_delta, max_treedepth, stepsize),
    seed    = seed
  )
}

# Prevalence recovery and diagnostics of an M6 fit; prefix is "m6d" or "m6j"
summarise_m6_sim8 = function(fit, dat_raw, prefix) {

  true_pr = dat_raw$mu_prevalence[1]
  draws = as.data.frame(fit)
  prev  = plogis(draws$b_pi_Intercept)

  out = tibble::tibble(
    prev_mean     = mean(prev),
    prev_lb       = quantile(prev, 0.025),
    prev_ub       = quantile(prev, 0.975),
    prev_bias     = mean(prev) - true_pr,
    prev_coverage = quantile(prev, 0.025) < true_pr &
      quantile(prev, 0.975) > true_pr,
    maxrhat       = max(rhat(fit), na.rm = TRUE),
    ndiv          = n_divergent(fit)
  )
  names(out) = paste0(prefix, "_", names(out))
  out
}

# Column order of the stored per-replicate results
sim8_result_columns = c(
  "sn", "mean_se", "mean_sp", "mu_prevalence", "tau_logit", "prop_gold", "kappa_prior",
  "disc_joint_map_se", "disc_joint_map_sp", "disc_marginal_se", "disc_marginal_sp",
  "disc_weighted_se", "disc_weighted_sp", "disc_top_weight",
  "disc_se_bias_wt", "disc_sp_bias_wt", "disc_sp_bias_marg",
  "kwga_prev_mean", "kwga_prev_lb", "kwga_prev_ub", "kwga_prev_bias",
  "kwga_prev_coverage", "kwga_n_eff_models",
  "joint_se_mean", "joint_sp_mean", "joint_se_sd", "joint_sp_sd",
  "joint_se_bias", "joint_sp_bias", "joint_se_coverage", "joint_sp_coverage",
  "joint_prev_mean", "joint_prev_bias", "joint_prev_coverage",
  "m6d_prev_mean", "m6d_prev_lb", "m6d_prev_ub", "m6d_prev_bias", "m6d_prev_coverage", "m6d_maxrhat",
  "m6j_prev_mean", "m6j_prev_lb", "m6j_prev_ub", "m6j_prev_bias", "m6j_prev_coverage", "m6j_maxrhat",
  "comp_maxrhat", "joint_maxrhat", "comp_ndiv", "joint_ndiv", "m6d_ndiv", "m6j_ndiv"
)

fit_sim8 = function(dat_raw,
                    se_grid = seq(0.60, 0.95, by = 0.025),
                    sp_grid = seq(0.60, 0.95, by = 0.025),
                    iter = 6000,
                    adapt_delta = 0.999,
                    max_treedepth = 20,
                    stepsize = NULL,
                    seed = NA) {

  # 1-2. Comparison model and discrete KWGA
  fit_comp = fit_comparison_sim8(dat_raw, iter = iter, adapt_delta = adapt_delta,
                                 max_treedepth = max_treedepth, stepsize = stepsize,
                                 seed = sim8_seed(seed, 1))
  kwga = estimate_sesp_kwga_sim8(fit_comp, dat_raw, se_grid = se_grid, sp_grid = sp_grid,
                                 seed = sim8_seed(seed, 2))

  # 3. Joint model
  fit_joint = fit_joint_sim8(dat_raw, iter = iter, adapt_delta = adapt_delta,
                             max_treedepth = max_treedepth, stepsize = stepsize,
                             seed = sim8_seed(seed, 3))
  joint = summarise_joint_sim8(fit_joint, dat_raw)

  # 4-5. M6 with KWGA-estimated (M6d) and joint-estimated (M6j) priors
  fit_m6d = fit_m6_from_estimates_sim8(dat_raw,
                                       se_est = kwga$disc_weighted_se, sp_est = kwga$disc_weighted_sp,
                                       iter = iter, adapt_delta = adapt_delta,
                                       max_treedepth = max_treedepth, stepsize = stepsize,
                                       seed = sim8_seed(seed, 4))
  fit_m6j = fit_m6_from_estimates_sim8(dat_raw,
                                       se_est = joint$joint_se_mean, sp_est = joint$joint_sp_mean,
                                       iter = iter, adapt_delta = adapt_delta,
                                       max_treedepth = max_treedepth, stepsize = stepsize,
                                       seed = sim8_seed(seed, 5))

  # 6. Collect metrics
  dplyr::bind_cols(
    tibble::tibble(
      sn            = dat_raw$sn[1],
      mean_se       = dat_raw$mean_se[1],
      mean_sp       = dat_raw$mean_sp[1],
      mu_prevalence = dat_raw$mu_prevalence[1],
      tau_logit     = dat_raw$tau_logit[1],
      prop_gold     = dat_raw$prop_gold[1],
      kappa_prior   = dat_raw$kappa_prior[1]
    ),
    kwga,
    joint,
    summarise_m6_sim8(fit_m6d, dat_raw, "m6d"),
    summarise_m6_sim8(fit_m6j, dat_raw, "m6j"),
    summarise_comparison_sim8(fit_comp)
  ) %>%
    dplyr::select(dplyr::all_of(sim8_result_columns))
}


# ---- Helper: Build M6 priors from estimated Se/Sp for bounded models -------------------------

# Build M6 priors using estimated (rather than external) Se/Sp centres.
# Gold measures get constant(10); screening measures get priors centred
# on the estimated values. Same as the prior function, except
# only for bounded models.

build_m6_priors_from_estimates_bounded = function(dat_raw,
                                                  priors_common,
                                                  se_est,
                                                  sp_est,
                                                  kappa        = 100,
                                                  floor_prob   = 0.5,
                                                  ceiling_prob = 0.999) {

  pri = c(priors_common$pi, priors_common$se_sp_sd)

  centers = dat_raw %>%
    dplyr::select(measure_id, is_gold) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  # Clamp external estimates (outer Se/Sp probability scale)
  se_c = clamp_probability(se_est, floor_prob, ceiling_prob)
  sp_c = clamp_probability(sp_est, floor_prob, ceiling_prob)

  for (i in seq_len(nrow(centers))) {
    m = centers$measure_id[i]

    if (centers$is_gold[i]) {
      pri = c(
        pri,
        set_prior("constant(10)", nlpar = "Se", class = "b",
                  coef = sprintf("measure_id%d", m)),
        set_prior("constant(10)", nlpar = "Sp", class = "b",
                  coef = sprintf("measure_id%d", m))
      )
    } else {
      pri = c(
        pri,
        prior_logit_for_measure(se_c, m, nlpar = "Se", n = kappa, bounded = TRUE),
        prior_logit_for_measure(sp_c, m, nlpar = "Sp", n = kappa, bounded = TRUE)
      )
    }
  }

  pri
}


# Helper: Get Divergent Transitions ---------------------------------------

# Total post-warmup divergent transitions for a brms fit: rstan stanfit path
# first, cmdstanr CmdStanMCMC fallback.
n_divergent <- function(fit) {
  tryCatch(
    sum(do.call(rbind, rstan::get_sampler_params(fit$fit, inc_warmup = FALSE))[, "divergent__"]),
    error = function(e) tryCatch(
      sum(fit$fit$sampler_diagnostics()[, , "divergent__"]),
      error = function(e2) NA_integer_)
  )
}

# Per-replicate seed for reproducible runs (NA = unseeded); cycle = 0 for the
# initial run and 1, 2, ... for successive repair cycles
sim8_replicate_seed = function(seed_base, s, r, cycle = 0) {
  if (is.na(seed_base)) NA else as.integer(seed_base + cycle * 1e8 + (s - 1) * 1e6 + r)
}

# Settings under which a stored result was fitted. Attached to every saved
# result as attr(x, "fit_settings") so that the archive records, per replicate,
# the seed, repair cycle and sampler settings that produced the file it holds.
sim8_fit_settings = function(seed, cycle, iter, adapt_delta, step_size = NULL) {
  list(seed = seed, cycle = cycle, iter = iter, adapt_delta = adapt_delta,
       step_size = if (is.null(step_size)) NA_real_ else step_size)
}

# ---- Runner: Loop over scenarios and replications -------------------------

run_sim8 = function(sim_data,
                    dir = "simulation_data/sim8/",
                    se_grid = seq(0.60, 0.95, by = 0.025),
                    sp_grid = seq(0.60, 0.95, by = 0.025),
                    iter = 6000,
                    adapt_delta = 0.999,
                    max_treedepth = 20,
                    skip_num = 0,
                    seed_base = NA) {

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  n_scen = length(sim_data)

  for (s in (seq_len(n_scen - skip_num) + skip_num)) {
    scen = sim_data[[s]]
    R = length(scen$replicates)

    for (r in seq_len(R)) {
      fn = file.path(dir, sprintf("Sim%d_%d.rds", s, r))

      if (file.exists(fn)) {
        message(sprintf("Skipping Sim %d, rep %d (exists)", s, r))
        next
      }

      message(sprintf("Sim %d/%d, Rep %d/%d", s, n_scen, r, R))

      dat_raw = scen$replicates[[r]]

      result = tryCatch({
        fit_sim8(
          dat_raw,
          se_grid = se_grid,
          sp_grid = sp_grid,
          iter = iter,
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth,
          seed = sim8_replicate_seed(seed_base, s, r)
        )
      }, error = function(e) {
        warning(sprintf("Error in scen %d rep %d: %s", s, r, e$message))
        NULL
      })

      if (!is.null(result)) {
        result$scenario_id = s
        result$rep_id = r
        attr(result, "fit_settings") = sim8_fit_settings(
          sim8_replicate_seed(seed_base, s, r), cycle = 0, iter = iter, adapt_delta = adapt_delta)
        saveRDS(result, fn)
      }
    }
  }

  message("Done.")
}


# ============================================================================
# Simulation 8: Repair Functions
# ============================================================================


# ---- Identify bad runs from saved results ---------------------------------

identify_bad_sim8 = function(dir = "simulation_data/sim8/") {
  files = list.files(dir, pattern = "Sim\\d+_\\d+\\.rds$", full.names = TRUE)

  if (length(files) == 0) {
    message("No result files found.")
    return(tibble::tibble(scenario_id = integer(), rep_id = integer(),
                          file = character()))
  }

  bad = purrr::map_dfr(files, function(fn) {
    res = tryCatch(readRDS(fn), error = function(e) NULL)
    if (is.null(res)) {
      # Corrupted file - treat as bad
      parts = regmatches(basename(fn),
                          regexec("Sim(\\d+)_(\\d+)", basename(fn)))[[1]]
      return(tibble::tibble(
        scenario_id = as.integer(parts[2]),
        rep_id = as.integer(parts[3]),
        file = fn,
        reason = "read_error"
      ))
    }

    # Check all model diagnostics
    is_bad = any(c(
      isTRUE(res$comp_maxrhat > 1.01),
      isTRUE(res$joint_maxrhat > 1.01),
      isTRUE(res$m6d_maxrhat > 1.01),
      isTRUE(res$m6j_maxrhat > 1.01),
      isTRUE(res$comp_ndiv  > 0),
      isTRUE(res$joint_ndiv > 0),
      isTRUE(res$m6d_ndiv   > 0),
      isTRUE(res$m6j_ndiv   > 0)
    ), na.rm = FALSE)

    if (is_bad) {
      tibble::tibble(
        scenario_id = res$scenario_id,
        rep_id = res$rep_id,
        file = fn,
        reason = dplyr::case_when(
          any(c(
            isTRUE(res$comp_maxrhat  > 1.01),
            isTRUE(res$joint_maxrhat > 1.01),
            isTRUE(res$m6d_maxrhat   > 1.01),
            isTRUE(res$m6j_maxrhat   > 1.01)
          )) ~ "rhat",

          any(c(
            isTRUE(res$comp_ndiv  > 0),
            isTRUE(res$joint_ndiv > 0),
            isTRUE(res$m6d_ndiv   > 0),
            isTRUE(res$m6j_ndiv   > 0)
          )) ~ "divergence",
          TRUE ~ "unknown"
        )
      )
    } else {
      NULL
    }
  })

  bad
}


# ---- Fix a set of bad runs -----------------------------------------------

fix_bad_sim8 = function(bad_trials,
                         sim_data,
                         dir = "simulation_data/sim8/",
                         cor_dir = "simulation_data/bad_runs/sim8/",
                         se_grid = seq(0.60, 0.95, by = 0.025),
                         sp_grid = seq(0.60, 0.95, by = 0.025),
                         iter = 12000,
                         adapt_delta = 0.999,
                         max_treedepth = 20,
                         step_size = 0.002,
                         seed_base = NA,
                         cycle = 1) {

  # Keep the previous result before replacing it with a repaired fit.
  if (!dir.exists(cor_dir) && !dir.create(cor_dir, recursive = TRUE)) {
    stop("Could not create backup directory: ", cor_dir)
  }

  nbad = nrow(bad_trials)

  for (i in seq_len(nbad)) {
    s = bad_trials$scenario_id[i]
    r = bad_trials$rep_id[i]

    fn = file.path(dir, sprintf("Sim%d_%d.rds", s, r))
    bad_fn = file.path(cor_dir, sprintf("Sim%d_%d.rds", s, r))

    message(sprintf("Fixing s=%d, r=%d  (%d of %d)", s, r, i, nbad))

    # Skip if already repaired in this cycle
    if (file.exists(bad_fn)) {
      message("    Already in corrected dir. Skipping.")
      next
    }

    # Copy the current result to the backup dir; it stays in place until the
    # repaired fit has been saved, so an interrupted or failed refit loses
    # nothing and the replicate remains visible to the next repair cycle
    if (file.exists(fn) && !file.copy(fn, bad_fn, overwrite = FALSE)) {
      stop("Could not back up ", fn, "; the original result has not been replaced.")
    }

    dat_raw = sim_data[[s]]$replicates[[r]]
    seed = sim8_replicate_seed(seed_base, s, r, cycle)
    message(sprintf("    repair cycle %d, seed %s", cycle, format(seed)))

    start_time = Sys.time()

    result = tryCatch({
      fit_sim8(
        dat_raw,
        se_grid = se_grid,
        sp_grid = sp_grid,
        iter = iter,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        stepsize = step_size,
        seed = seed
      )
    }, error = function(e) {
      warning(sprintf("Error repairing s=%d r=%d: %s", s, r, e$message))
      NULL
    })

    elapsed = as.numeric(Sys.time() - start_time, units = "secs")

    if (!is.null(result)) {
      result$scenario_id = s
      result$rep_id = r
      attr(result, "fit_settings") = sim8_fit_settings(
        seed, cycle = cycle, iter = iter, adapt_delta = adapt_delta, step_size = step_size)
      tmp_fn = paste0(fn, ".tmp")
      saveRDS(result, tmp_fn)
      if (!file.rename(tmp_fn, fn)) stop("Could not replace ", fn, " with the repaired result.")
      message(sprintf("    Completed s=%d, r=%d in %.1f s", s, r, elapsed))
    } else {
      message(sprintf("    FAILED s=%d, r=%d after %.1f s (previous result kept)", s, r, elapsed))
    }
  }
}


# ---- Repair cycle ---------------------------------------------------------

repair_cycle_sim8 = function(sim_data,
                             dir = "simulation_data/sim8/",
                             cor_dir = "simulation_data/bad_runs/sim8/",
                             se_grid = seq(0.60, 0.95, by = 0.025),
                             sp_grid = seq(0.60, 0.95, by = 0.025),
                             iter = 12000,
                             adapt_delta = 0.999,
                             max_treedepth = 20,
                             step_size = 0.002,
                             cycles = 3,
                             scenario_ids = NULL,
                             seed_base = NA) {

  for (cyc in seq_len(cycles)) {
    message(sprintf("\n=== Repair cycle %d / %d ===", cyc, cycles))

    bad = identify_bad_sim8(dir = dir)

    # Filter to assigned scenarios if partitioned
    if (!is.null(scenario_ids)) {
      bad = bad %>% dplyr::filter(scenario_id %in% scenario_ids)
    }

    if (nrow(bad) == 0) {
      message("No bad runs found in assigned scenarios. Done.")
      return(0L)
    }

    message(sprintf("Found %d bad runs", nrow(bad)))

    for (j in seq_len(nrow(bad))) {
      bad_fn = file.path(cor_dir,
                         sprintf("Sim%d_%d.rds",
                                 bad$scenario_id[j], bad$rep_id[j]))
      if (file.exists(bad_fn)) file.remove(bad_fn)
    }

    fix_bad_sim8(
      bad_trials = bad,
      sim_data = sim_data,
      dir = dir,
      cor_dir = cor_dir,
      se_grid = se_grid,
      sp_grid = sp_grid,
      iter = iter,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      step_size = step_size,
      seed_base = seed_base,
      cycle = cyc
    )
  }

  bad_final = identify_bad_sim8(dir = dir)
  if (!is.null(scenario_ids)) {
    bad_final = bad_final %>% dplyr::filter(scenario_id %in% scenario_ids)
  }
  n_remaining = nrow(bad_final)

  if (n_remaining > 0) {
    message(sprintf("\n%d bad runs remain after %d cycles:", n_remaining, cycles))
    print(bad_final %>% dplyr::select(scenario_id, rep_id, reason))
  } else {
    message("\nAll assigned runs clean.")
  }

  n_remaining
}


# ---- Read results from per-replicate files --------------------------------

read_sim8_rds = function(dir = "simulation_data/sim8/") {
  files = list.files(dir, pattern = "Sim\\d+_\\d+\\.rds$", full.names = TRUE)

  if (length(files) == 0) {
    warning("No result files found in ", dir)
    return(NULL)
  }

  purrr::map_dfr(files, function(fn) {
    tryCatch(readRDS(fn), error = function(e) {
      warning(sprintf("Could not read %s: %s", fn, e$message))
      NULL
    })
  })
}


# ---- Process Simulation 8 ------------------------------------------------

process_sim8 = function(
    dir,
    group_vars,
    exclude_bad = TRUE
) {

  results = read_sim8_rds(dir = dir)

  if (is.null(results) || nrow(results) == 0) {
    stop("No results found in ", dir)
  }

  # Identify bad replications.
  # A rep is "bad" if ANY sub-model has Rhat > 1.01 OR any divergence;
  # coalesce(., TRUE) treats unverifiable diagnostics (NA) as bad.
  results = results %>%
    dplyr::mutate(.is_bad = dplyr::coalesce(
      comp_maxrhat > 1.01 | joint_maxrhat > 1.01 |
        m6d_maxrhat > 1.01 | m6j_maxrhat > 1.01 |
        comp_ndiv > 0 | joint_ndiv > 0 |
        m6d_ndiv > 0 | m6j_ndiv > 0,
      TRUE))

  bad_reps = results %>%
    dplyr::filter(.is_bad) %>%
    dplyr::select(scenario_id, rep_id) %>%
    dplyr::distinct()

  # Per-cell exclusion rate + 5% cap (per FULL design cell = scenario_id),
  # computed BEFORE exclusion so it reflects what is actually dropped.
  cell_exclusion = results %>%
    dplyr::group_by(scenario_id) %>%
    dplyr::summarize(n_total  = dplyr::n(),
                     n_bad    = sum(.is_bad),
                     prop_bad = mean(.is_bad),
                     .groups  = "drop") %>%
    dplyr::mutate(over_cap = prop_bad > 0.05)

  if (any(cell_exclusion$over_cap)) {
    warning(sprintf(
      "exclude_bad: >5%% of reps flagged in %d/%d cells (scenarios %s).",
      sum(cell_exclusion$over_cap), nrow(cell_exclusion),
      paste(cell_exclusion$scenario_id[cell_exclusion$over_cap],
            collapse = ", ")))
  }

  if (exclude_bad) {
    results = results %>% dplyr::filter(!.is_bad)
  }
  results = dplyr::select(results, -.is_bad)

  # Summarize
  results_grouped = results %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars)))

  summary = results_grouped %>%
    dplyr::summarize(

      # Discrete KWGA: Sp recovery
      disc_sp_mae_wt      = mean(abs(disc_sp_bias_wt), na.rm = TRUE),
      disc_sp_mae_marg    = mean(abs(disc_sp_bias_marg), na.rm = TRUE),
      disc_sp_bias_wt     = mean(disc_sp_bias_wt, na.rm = TRUE),
      disc_sp_bias_marg   = mean(disc_sp_bias_marg, na.rm = TRUE),

      # Discrete KWGA: Se recovery
      disc_se_mae_wt      = mean(abs(disc_se_bias_wt), na.rm = TRUE),
      disc_se_bias_wt     = mean(disc_se_bias_wt, na.rm = TRUE),

      # Discrete KWGA: concentration
      disc_top_weight     = mean(disc_top_weight, na.rm = TRUE),
      disc_n_eff_models   = mean(kwga_n_eff_models, na.rm = TRUE),

      # KWGA prevalence
      kwga_prev_rmse      = sqrt(mean(kwga_prev_bias^2, na.rm = TRUE)),
      kwga_prev_bias      = mean(kwga_prev_bias, na.rm = TRUE),
      kwga_prev_coverage  = mean(kwga_prev_coverage, na.rm = TRUE),

      # Joint model: Sp recovery
      joint_sp_mae        = mean(abs(joint_sp_bias), na.rm = TRUE),
      joint_sp_bias       = mean(joint_sp_bias, na.rm = TRUE),
      joint_sp_coverage   = mean(joint_sp_coverage, na.rm = TRUE),
      joint_sp_sd         = mean(joint_sp_sd, na.rm = TRUE),

      # Joint model: Se recovery
      joint_se_mae        = mean(abs(joint_se_bias), na.rm = TRUE),
      joint_se_bias       = mean(joint_se_bias, na.rm = TRUE),
      joint_se_coverage   = mean(joint_se_coverage, na.rm = TRUE),
      joint_se_sd         = mean(joint_se_sd, na.rm = TRUE),

      # Joint model: prevalence
      joint_prev_rmse     = sqrt(mean(joint_prev_bias^2, na.rm = TRUE)),
      joint_prev_bias     = mean(joint_prev_bias, na.rm = TRUE),
      joint_prev_coverage = mean(joint_prev_coverage, na.rm = TRUE),

      # M6 with discrete priors: prevalence
      m6d_prev_rmse       = sqrt(mean(m6d_prev_bias^2, na.rm = TRUE)),
      m6d_prev_bias       = mean(m6d_prev_bias, na.rm = TRUE),
      m6d_prev_coverage   = mean(m6d_prev_coverage, na.rm = TRUE),
      m6d_ci_width        = mean(m6d_prev_ub - m6d_prev_lb, na.rm = TRUE),

      # M6 with joint priors: prevalence
      m6j_prev_rmse       = sqrt(mean(m6j_prev_bias^2, na.rm = TRUE)),
      m6j_prev_bias       = mean(m6j_prev_bias, na.rm = TRUE),
      m6j_prev_coverage   = mean(m6j_prev_coverage, na.rm = TRUE),
      m6j_ci_width        = mean(m6j_prev_ub - m6j_prev_lb, na.rm = TRUE),

      # Diagnostics
      prop_bad_comp       = mean(comp_maxrhat > 1.01, na.rm = TRUE),
      prop_bad_joint      = mean(joint_maxrhat > 1.01, na.rm = TRUE),
      prop_bad_m6d        = mean(m6d_maxrhat > 1.01, na.rm = TRUE),
      prop_bad_m6j        = mean(m6j_maxrhat > 1.01, na.rm = TRUE),

      n = dplyr::n(),
      .groups = "drop"
    )

  list(
    summary        = summary,
    bad_reps       = bad_reps,
    n_bad          = nrow(bad_reps),
    cell_exclusion = cell_exclusion
  )
}
