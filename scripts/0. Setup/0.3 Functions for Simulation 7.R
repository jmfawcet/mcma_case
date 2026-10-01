# ============================================================================
# Simulation 7: Moderator Power
# ============================================================================
#
# Evaluates whether including corrected screening studies improves the
# detection of moderator effects compared to gold-standard-only analysis.
#
# Design: 3 (moderator slope: 0, 0.5, 1.0) x
#         3 (proportion gold: .10, .25, .50) x
#         2 (prior centre: .85/.75, .85/.85) = 18 scenarios
#
# Six models per replication:
#   M1:        Naive meta-regression (ignores misclassification)
#   Gold-only: Meta-regression on gold standard studies only
#   KnownSeSp: All studies corrected deterministically with their true Se/Sp
#   AllGold:   Every study re-observed as y* ~ Binomial(n, theta_true)
#   M6:        Corrected meta-regression (screening priors, gold fixed)
#   M8:        Joint model estimating screening Se/Sp
#
# ============================================================================


# ---- Scenario Grid --------------------------------------------------------

scenario_grid_sim7 = function(
    sn_levels       = 40,
    accuracy        = 0.85,      # true mu_se
    sp_offset       = 0.10,      # true mu_sp = accuracy - sp_offset
    sigma_levels    = 0.2,
    omega_levels    = 0.2,
    mu_prev_vals    = 0.03,      # baseline prevalence
    tau             = 1.0,
    rho             = -0.5,
    kappa_prior     = 100,
    slope_vals      = c(0, 0.5, 1.0),
    prop_gold_vals  = c(0.10, 0.25, 0.50),
    prior_se_vals   = c(0.85, 0.85),
    prior_sp_vals   = c(0.75, 0.85)
) {
  # Build prior configurations as paired vectors
  prior_configs = data.frame(prior_se = prior_se_vals, prior_sp = prior_sp_vals)

  base_grid = expand.grid(
    sn          = sn_levels,
    acc_level   = accuracy,
    sp_offset   = sp_offset,
    sig_level   = sigma_levels,
    omg_level   = omega_levels,
    mu_prev     = mu_prev_vals,
    tau         = tau,
    rho         = rho,
    kappa_prior = kappa_prior,
    slope       = slope_vals,
    prop_gold   = prop_gold_vals,
    prior_idx   = seq_len(nrow(prior_configs))
  )

  base_grid$prior_se = prior_configs$prior_se[base_grid$prior_idx]
  base_grid$prior_sp = prior_configs$prior_sp[base_grid$prior_idx]
  base_grid$prior_idx = NULL

  base_grid
}


# ---- Data Generator with Moderator and Gold Standard ----------------------

simulate_dataset_moderator = function(
    sn = 40,
    meanlog_n = 5.50, sdlog_n = 0.55, add_const_n = 30,
    di = 10,
    mean_se = 0.85, mean_sp = 0.75,
    sd_se_logit = 0.2, sd_sp_logit = 0.2, rho = -0.5,
    omega_se = 0.2, omega_sp = 0.2,
    mu_prevalence = 0.03,
    tau_logit = 1.0,
    kappa_prior = 100,
    floor_prob = 0.5, ceiling_prob = 0.999,
    prop_gold = 0.10,
    slope = 0,
    prior_se = 0.85,
    prior_sp = 0.75
) {
  # Determine gold vs screening
  n_gold   = max(1, round(sn * prop_gold))
  n_screen = sn - n_gold

  # Sample sizes for all studies
  n_all = sample_sizes_logn(sn, meanlog_n, sdlog_n, add_const_n)

  # Continuous moderator ~ N(0, 1)
  x_mod = rnorm(sn, mean = 0, sd = 1)

  # True prevalence with moderator effect:
  # logit(theta_i) = beta0 + beta1 * x_i + u_i
  beta0 = qlogis(mu_prevalence)
  u_i   = rnorm(sn, mean = 0, sd = tau_logit)
  theta_all = plogis(beta0 + slope * x_mod + u_i)

  # Gold standard studies (indices 1:n_gold)
  Se_gold = rep(ceiling_prob, n_gold)
  Sp_gold = rep(ceiling_prob, n_gold)

  # Screening studies (indices (n_gold+1):sn)
  measures = draw_measures(di, mean_se, mean_sp, sd_se_logit, sd_sp_logit,
                           rho, floor_prob, ceiling_prob)
  measure_idx = sample(measures$measure_id, size = n_screen, replace = TRUE)

  acc = jitter_study_accuracy(measures, measure_idx,
                              omega_se = omega_se, omega_sp = omega_sp,
                              floor_prob = floor_prob, ceiling_prob = ceiling_prob)

  # Combine Se/Sp
  Se_all = c(Se_gold, acc$Se)
  Sp_all = c(Sp_gold, acc$Sp)

  # Forward misclassification for all studies
  fm = forward_misclassification(n_all, theta_all, Se_all, Sp_all)

  # Measure IDs: 0 for gold, actual IDs for screening
  measure_id_all = c(rep(0, n_gold), measure_idx)

  # Base Se/Sp for priors
  se_base_all = c(rep(ceiling_prob, n_gold), measures$se_base[measure_idx])
  sp_base_all = c(rep(ceiling_prob, n_gold), measures$sp_base[measure_idx])
  eta_se_base_all = c(rep(qlogis(ceiling_prob), n_gold), measures$eta_se_base[measure_idx])
  eta_sp_base_all = c(rep(qlogis(ceiling_prob), n_gold), measures$eta_sp_base[measure_idx])

  data.frame(
    sn            = sn,
    di            = di,
    mean_se       = mean_se,
    mean_sp       = mean_sp,
    sd_se_logit   = sd_se_logit,
    sd_sp_logit   = sd_sp_logit,
    rho           = rho,
    omega_se      = omega_se,
    omega_sp      = omega_sp,
    sigma_bias    = 0,
    kappa_prior   = kappa_prior,
    mu_prevalence = mu_prevalence,
    tau_logit     = tau_logit,
    slope_true    = slope,
    prop_gold     = prop_gold,
    prior_se      = prior_se,
    prior_sp      = prior_sp,
    study         = seq_len(sn),
    n             = n_all,
    y             = fm$y,
    p_obs         = fm$p_obs,
    theta_true    = theta_all,
    x_mod         = x_mod,
    Se            = Se_all,
    Sp            = Sp_all,
    measure_id    = measure_id_all,
    se_base       = se_base_all,
    sp_base       = sp_base_all,
    eta_se_base   = eta_se_base_all,
    eta_sp_base   = eta_sp_base_all,
    is_gold       = c(rep(TRUE, n_gold), rep(FALSE, n_screen)),
    stringsAsFactors = FALSE # Fixed an issue I was having, not sure if it is needed though
  )
}


# ---- Simulation Generator -------------------------------------------------

simulate_all_scenarios_sim7 = function(R = 200,
                                       base_seed = 20250903,
                                       meanlog_n = 5.50, sdlog_n = 0.55,
                                       add_const_n = 30,
                                       di = 10,
                                       floor_prob = 0.5,
                                       ceiling_prob = 0.999,
                                       grid = scenario_grid_sim7()) {
  n_scen = nrow(grid)
  out = vector("list", length = n_scen)

  for (s in seq_len(n_scen)) {
    message(sprintf("Scenario %d / %d", s, n_scen))

    g = grid[s, ]

    mu_se = g$acc_level
    mu_sp = g$acc_level - g$sp_offset

    datasets = vector("list", length = R)
    for (r in seq_len(R)) {
      set.seed(base_seed + (s - 1) * 1e6 + r)
      datasets[[r]] = simulate_dataset_moderator(
        sn = g$sn,
        meanlog_n = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
        di = di,
        mean_se = mu_se, mean_sp = mu_sp,
        sd_se_logit = g$sig_level, sd_sp_logit = g$sig_level, rho = g$rho,
        omega_se = g$omg_level, omega_sp = g$omg_level,
        mu_prevalence = g$mu_prev,
        tau_logit = g$tau,
        kappa_prior = g$kappa_prior,
        floor_prob = floor_prob, ceiling_prob = ceiling_prob,
        prop_gold = g$prop_gold,
        slope = g$slope,
        prior_se = g$prior_se,
        prior_sp = g$prior_sp
      )
    }

    out[[s]] = list(
      scenario_id = s,
      factors = list(
        sn = g$sn, accuracy = mu_se, mu_prev = g$mu_prev,
        rho = g$rho, slope = g$slope, prop_gold = g$prop_gold,
        prior_se = g$prior_se, prior_sp = g$prior_sp,
        kappa_prior = g$kappa_prior
      ),
      replicates = datasets
    )
  }
  out
}


# ---- Prior Builder: Gold + Common Centre ----------------------------------

# Build M6 priors for moderator model with gold + screening.
# Gold measures (measure_id == 0): constant(10).
# Screening measures: single common prior centre, no sigma_bias.
# This is simplified and uses some defaults specific to this simulation,
# but basically the same as before for M6.

build_priors_m6_mod = function(dat_fac,
                               priors_common,
                               n_eff        = 100,
                               prior_se     = 0.85,
                               prior_sp     = 0.75,
                               floor_prob   = 0.5,
                               ceiling_prob = 0.999) {

  pri = c(priors_common$pi_mod, priors_common$se_sp_sd)

  centers = dat_fac %>%
    dplyr::select(measure_id, is_gold) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  gold_centers   = centers %>% dplyr::filter(is_gold)
  screen_centers = centers %>% dplyr::filter(!is_gold)

  # Gold: fixed near-perfect
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

  # Screening: common centre, no noise
  se_common = clamp_probability(prior_se, floor_prob, ceiling_prob)
  sp_common = clamp_probability(prior_sp, floor_prob, ceiling_prob)

  for (i in seq_len(nrow(screen_centers))) {
    m = screen_centers$measure_id[i]
    pri = c(
      pri,
      prior_logit_for_measure(se_common, m, nlpar = "Se", n = n_eff),
      prior_logit_for_measure(sp_common, m, nlpar = "Sp", n = n_eff)
    )
  }

  pri
}


# ---- Common Priors with Moderator -----------------------------------------

# Extends make_common_priors to include the moderator slope prior
make_common_priors_mod = function(true_pr) {
  list(
    # M1 / Gold-only / AllGold: intercept + moderator + random intercept
    m1_mod = c(
      set_prior(sprintf("normal(%f, 1.5)", true_pr), class = "Intercept"),
      prior(normal(0, 1), class = "b", coef = "x_mod"),
      prior(normal(0, 1), class = "sd")
    ),
    # M6: pi submodel with intercept + moderator + random intercept
    pi_mod = c(
      set_prior(sprintf("normal(%f, 1.5)", true_pr), nlpar = "pi"),
      prior(normal(0, 1), nlpar = "pi", class = "b", coef = "x_mod"),
      prior(normal(0, 1), nlpar = "pi", class = "sd")
    ),
    # SD priors for Se/Sp study REs
    se_sp_sd = c(
      prior(normal(0, .5), nlpar = "Se", class = "sd"),
      prior(normal(0, .5), nlpar = "Sp", class = "sd")
    )
  )
}


# ---- Prior Builder: M8 with moderator ------------------------------------

build_priors_m8_mod <- function(true_pr,
                                se_prior_center = 0.80,
                                sp_prior_center = 0.80,
                                se_prior_sd     = 1.0,
                                sp_prior_sd     = 1.0) {

  se_inner <- qlogis(map_to_inner_t(se_prior_center))
  sp_inner <- qlogis(map_to_inner_t(sp_prior_center))

  c(
    # Prevalence intercept + moderator + RE
    set_prior(sprintf("normal(%f, 1.5)", qlogis(true_pr)), nlpar = "pi"),
    prior(normal(0, 1), nlpar = "pi", class = "b", coef = "x_mod"),
    prior(normal(0, 1), nlpar = "pi", class = "sd"),

    # Se/Sp RE SDs
    prior(normal(0, 0.5), nlpar = "Se", class = "sd"),
    prior(normal(0, 0.5), nlpar = "Sp", class = "sd"),

    # Gold: fixed near-perfect
    prior_string("constant(10)", class = "b",
                 coef = "acc_groupgold", nlpar = "Se"),
    prior_string("constant(10)", class = "b",
                 coef = "acc_groupgold", nlpar = "Sp"),

    # Screening: wide priors on inner logit scale
    prior_string(sprintf("normal(%f, %f)", se_inner, se_prior_sd),
                 class = "b", coef = "acc_groupscreen", nlpar = "Se"),
    prior_string(sprintf("normal(%f, %f)", sp_inner, sp_prior_sd),
                 class = "b", coef = "acc_groupscreen", nlpar = "Sp")
  )
}


# ---- brms Formulas with Moderator -----------------------------------------

bf_m6_mod = bf(
  y | trials(n) ~ (inv_logit(pi) * inv_logit(Se)) +
    ((1 - inv_logit(pi)) * (1 - inv_logit(Sp))),
  pi ~ 1 + x_mod + (1 | s | study),
  Se ~ measure_id - 1 + (1 | s | study),
  Sp ~ measure_id - 1 + (1 | s | study),
  nl = TRUE
)

# ---- Formula: M8 with moderator (bounded) --------------------------------

bf_m8_mod <- bf(
  y | trials(n) ~ inv_logit(pi) * (0.5 + 0.5 * inv_logit(Se)) +
    (1 - inv_logit(pi)) * (1 - (0.5 + 0.5 * inv_logit(Sp))),
  pi ~ 1 + x_mod + (1 | s | study),
  Se ~ 0 + acc_group + (1 | s | study),
  Sp ~ 0 + acc_group + (1 | s | study),
  nl = TRUE
)


# ---- Moderator Metrics ----------------------------------------------------

# Extract moderator slope metrics from a brmsfit
extract_slope_metrics = function(fit, dat_raw, model_name, is_nlpar = FALSE) {

  true_slope = dat_raw$slope_true[1]
  dr = as_draws_df(fit)

  # Find the slope draws
  if (is_nlpar) {
    slope_draws = dr$b_pi_x_mod
  } else {
    slope_draws = dr$b_x_mod
  }

  if (is.null(slope_draws)) {
    warning(sprintf("Could not find slope draws for %s", model_name))
    return(NULL)
  }

  slope_mean = mean(slope_draws)
  slope_lb   = quantile(slope_draws, 0.025)
  slope_ub   = quantile(slope_draws, 0.975)

  # Prevalence intercept
  if (is_nlpar) {
    prev_draws = plogis(dr$b_pi_Intercept)
  } else {
    prev_draws = plogis(dr$b_Intercept)
  }

  prev_mean = mean(prev_draws)
  prev_lb   = quantile(prev_draws, 0.025)
  prev_ub   = quantile(prev_draws, 0.975)
  true_prev = dat_raw$mu_prevalence[1]

  tibble(
    model         = model_name,

    # Slope metrics
    slope_mean    = slope_mean,
    slope_sd      = sd(slope_draws),
    slope_lb      = slope_lb,
    slope_ub      = slope_ub,
    slope_bias    = slope_mean - true_slope,
    slope_sqerror = (slope_mean - true_slope)^2,
    slope_coverage = slope_lb < true_slope & slope_ub > true_slope,
    slope_power   = (slope_lb > 0) | (slope_ub < 0),  # CI excludes zero
    slope_width   = slope_ub - slope_lb,

    # Prevalence intercept metrics
    prev_mean     = prev_mean,
    prev_lb       = prev_lb,
    prev_ub       = prev_ub,
    prev_bias     = prev_mean - true_prev,
    prev_sqerror  = (prev_mean - true_prev)^2,
    prev_coverage = prev_lb < true_prev & prev_ub > true_prev,

    # Sampler diagnostics
    maxrhat       = max(rhat(fit), na.rm = TRUE),
    ndiv          = sum(do.call(rbind,
                                rstan::get_sampler_params(fit$fit, inc_warmup = FALSE))[, "divergent__"])
  )
}


# ---- Slope Metrics Extractor for M8 --------------------------------------

extract_slope_metrics_m8 <- function(fit, dat_raw) {

  true_slope <- dat_raw$slope_true[1]
  true_prev  <- dat_raw$mu_prevalence[1]
  dr <- as_draws_df(fit)

  # Slope

  slope_draws <- dr$b_pi_x_mod
  slope_mean  <- mean(slope_draws)
  slope_lb    <- quantile(slope_draws, 0.025)
  slope_ub    <- quantile(slope_draws, 0.975)

  # Prevalence intercept
  prev_draws <- plogis(dr$b_pi_Intercept)
  prev_mean  <- mean(prev_draws)
  prev_lb    <- quantile(prev_draws, 0.025)
  prev_ub    <- quantile(prev_draws, 0.975)

  # Se/Sp (bounded back-transform)
  se_draws <- 0.5 + 0.5 * plogis(dr$b_Se_acc_groupscreen)
  sp_draws <- 0.5 + 0.5 * plogis(dr$b_Sp_acc_groupscreen)

  tibble(
    model          = "M8",

    # Slope metrics (same as Sim 7 output)
    slope_mean     = slope_mean,
    slope_sd       = sd(slope_draws),
    slope_lb       = slope_lb,
    slope_ub       = slope_ub,
    slope_bias     = slope_mean - true_slope,
    slope_sqerror  = (slope_mean - true_slope)^2,
    slope_coverage = slope_lb < true_slope & slope_ub > true_slope,
    slope_power    = (slope_lb > 0) | (slope_ub < 0),
    slope_width    = slope_ub - slope_lb,

    # Prevalence intercept metrics
    prev_mean      = prev_mean,
    prev_lb        = prev_lb,
    prev_ub        = prev_ub,
    prev_bias      = prev_mean - true_prev,
    prev_sqerror   = (prev_mean - true_prev)^2,
    prev_coverage  = prev_lb < true_prev & prev_ub > true_prev,

    # Se/Sp recovery
    se_mean        = mean(se_draws),
    sp_mean        = mean(sp_draws),
    se_bias        = mean(se_draws) - dat_raw$mean_se[1],
    sp_bias        = mean(sp_draws) - dat_raw$mean_sp[1],

    # Diagnostics
    maxrhat        = max(rhat(fit), na.rm = TRUE),
    ndiv           = sum(do.call(rbind,
                                 rstan::get_sampler_params(fit$fit,
                                                           inc_warmup = FALSE))[, "divergent__"])
  )
}


# ---- Fitting Function -----------------------------------------------------

fit_models_sim7 = function(dat_raw,
                           iter = 6000,
                           adapt_delta = 0.999,
                           max_treedepth = 20,
                           stepsize = NULL) {

  true_pr    = dat_raw$mu_prevalence[1]
  kappa_pr   = dat_raw$kappa_prior[1]
  c_prior_se = dat_raw$prior_se[1]
  c_prior_sp = dat_raw$prior_sp[1]

  common_priors = make_common_priors_mod(qlogis(true_pr))

  if (is.null(stepsize)) {
    ctrl = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  } else {
    ctrl = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                step_size = stepsize)
  }

  # Normalize scenario metadata (gold rows may differ)
  scenario_cols = c("sn", "di", "mean_se", "mean_sp", "sd_se_logit",
                    "sd_sp_logit", "rho", "omega_se", "omega_sp")
  screen_row = dat_raw %>% dplyr::filter(!is_gold) %>% dplyr::slice(1)
  for (col in scenario_cols) {
    dat_raw[[col]] = screen_row[[col]]
  }

  # M1: Naive meta-regression (all studies, ignores misclassification)
  dat_m1 = dat_raw %>% dplyr::select(y, n, study, x_mod)

  fit1 = brm(
    y | trials(n) ~ 1 + x_mod + (1 | study),
    data    = dat_m1,
    family  = binomial(),
    backend = 'cmdstanr',
    prior   = common_priors$m1_mod,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # Gold-only: Meta-regression on gold standard studies only (precise number wil vary depending on p(gold))
  dat_gold = dat_raw %>% dplyr::filter(is_gold) %>% dplyr::select(y, n, study, x_mod)

  gold_priors = c(
    set_prior(sprintf("normal(%f, 1.5)", qlogis(true_pr)), class = "Intercept"),
    prior(normal(0, 1), class = "b", coef = "x_mod"),
    prior(normal(0, 1), class = "sd")
  )

  fit_gold = brm(
    y | trials(n) ~ 1 + x_mod + (1 | study),
    data    = dat_gold,
    family  = binomial(),
    backend = 'cmdstanr',
    prior   = gold_priors,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # Known Se/Sp: all studies corrected deterministically (Rogan-Gladen) using their true Se/Sp
  dat_known = dat_raw %>%
    dplyr::mutate(y = round(n * rogan_gladen(y / n, Se, Sp))) %>%
    dplyr::select(y, n, study, x_mod)

  fit_known = brm(
    y | trials(n) ~ 1 + x_mod + (1 | study),
    data    = dat_known,
    family  = binomial(),
    backend = 'cmdstanr',
    prior   = common_priors$m1_mod,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # AllGold: oracle benchmark. Every study is re-observed as y* ~ Binomial(n, theta_true)
  # and the naive meta-regression is fitted to y*. The draw is made at fitting time and
  # is not seeded, so AllGold rows vary slightly between runs.
  dat_allgold = dat_raw %>%
    dplyr::mutate(y = rbinom(dplyr::n(), size = n, prob = theta_true)) %>%
    dplyr::select(y, n, study, x_mod)

  fit_allgold = brm(
    y | trials(n) ~ 1 + x_mod + (1 | study),
    data    = dat_allgold,
    family  = binomial(),
    backend = 'cmdstanr',
    prior   = common_priors$m1_mod,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # M6: Corrected meta-regression with moderator
  dat_m6 = dat_raw %>%
    dplyr::select(y, n, study, x_mod, measure_id, se_base, sp_base, is_gold)

  pri6 = build_priors_m6_mod(
    dat_m6, common_priors,
    n_eff    = kappa_pr,
    prior_se = c_prior_se,
    prior_sp = c_prior_sp
  )

  fit6 = brm(
    formula = bf_m6_mod,
    data    = dat_m6 %>% dplyr::mutate(measure_id = as.character(measure_id)),
    family  = binomial(link = 'identity'),
    backend = 'cmdstanr',
    prior   = pri6,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # M8: Joint model estimating screening Se/Sp (gold studies fixed near-perfect)
  dat_m8 = dat_raw %>%
    dplyr::select(y, n, study, x_mod, is_gold) %>%
    dplyr::mutate(acc_group = factor(ifelse(is_gold, "gold", "screen")))
  pri8 = build_priors_m8_mod(true_pr)

  fit8 = brm(
    formula = bf_m8_mod,
    data    = dat_m8,
    family  = binomial(link = "identity"),
    backend = 'cmdstanr',
    prior   = pri8,
    chains  = 4, cores = 4,
    iter    = iter,
    refresh = 0, silent = 2,
    control = ctrl
  )

  # Collect metrics
  k = dplyr::bind_rows( # Still feeling k, it is tradition now
    extract_slope_metrics(fit1, dat_raw, "M1", is_nlpar = FALSE),
    extract_slope_metrics(fit_gold, dat_raw, "Gold", is_nlpar = FALSE),
    extract_slope_metrics(fit_known, dat_raw, "KnownSeSp", is_nlpar = FALSE),
    extract_slope_metrics(fit_allgold, dat_raw, "AllGold", is_nlpar = FALSE),
    extract_slope_metrics(fit6, dat_raw, "M6", is_nlpar = TRUE),
    extract_slope_metrics_m8(fit8, dat_raw)
  ) %>%
    dplyr::bind_cols(
      dat_raw %>%
        dplyr::select(sn, di, mean_se, mean_sp, mu_prevalence, tau_logit,
                      slope_true, prop_gold, prior_se, prior_sp, kappa_prior) %>%
        dplyr::distinct() %>%
        dplyr::slice(1)
    )

  return(k)
}


# ---- Run / Repair Functions -----------------------------------------------

run_all_rds_sim7 = function(sims, dir,
                            iter = 6000, adapt_delta = 0.999,
                            skip_num = 0) {
  for (i in (seq_len(length(sims) - skip_num) + skip_num)) {
    message(sprintf("Starting scenario i = %d", i))
    c_sim = sims[[i]]

    for (j in seq_len(length(c_sim$replicates))) {
      fn = sprintf('%sSim%d_%d.rds', dir, i, j)

      if (file.exists(fn)) {
        message(sprintf("  Loading r = %d", j))
        next
      }

      start_time = Sys.time()
      c_rep = c_sim$replicates[[j]]

      summ = fit_models_sim7(c_rep, iter = iter, adapt_delta = adapt_delta)

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

repair_cycle_sim7 = function(sims,
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

    for (j in seq_len(nrow(bad_trials))) {
      bad_fn = sprintf('%sSim%d_%d.rds', cor_dir,
                       bad_trials$sim_id[j], bad_trials$rep_id[j])
      if (file.exists(bad_fn)) file.remove(bad_fn)
    }

    if (nrow(bad_trials) > 0) {
      fix_all_rds_sim7(bad_trials, sims,
                       iter = iter, adapt_delta = adapt_delta,
                       step_size = step_size,
                       dir = dir, cor_dir = cor_dir)
    } else {
      return(0)
    }
  }

  results_dat = read_sim_rds(dir = dir)
  bad_trials = results_dat %>%
    dplyr::filter(ndiv > 0 | maxrhat > 1.01) %>%
    dplyr::select(sim_id, rep_id) %>%
    dplyr::distinct()
  return(nrow(bad_trials))
}

fix_all_rds_sim7 = function(bad_trials, sims,
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
    summ = fit_models_sim7(dat_raw, iter = iter, adapt_delta = adapt_delta,
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
