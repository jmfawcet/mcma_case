# ============================================================================
# General Libraries and Functions
# ============================================================================
#
# These were the first functions I created, which are used throughout my code.
# When I wrote them, I hadn't imagined some of the simulations yet, so I ended up
# with variants of some of these, but they are the 'core' of the simulations and
# case study. Some libraries may no longer be needed.
# ============================================================================

# Libraries ---------------------------------------------------------------

library(MASS)
library(patchwork)
library(gghalves)
library(posterior)
library(tibble)
library(rlang)
library(openxlsx)
library(cmdstanr)
library(brms)
library(purrr)
library(tidyr)
library(tidyverse)

# Simulation Functions ----------------------------------------------------

## Rogan-Gladen Correction Helper Function ---------------------------------

rogan_gladen = function(p_obs, Se, Sp, clamp = TRUE) {
  # Basic RG correction
  p_corr = (p_obs + Sp - 1) / (Se + Sp - 1)

  # If the correction goes above 1 or below 0, set it to 1 or 0
  if (clamp) {
    p_corr = clamp_probability(p_corr, 0, 1)
  }

  return(p_corr)
}

## Generate Measures -------------------------------------------------------

# Returns measure-level Se/Sp (in both logit and probability space)
draw_measures = function(di = 10,  # Number of measures
                          mean_se = 0.80, mean_sp = 0.75,  # Mean Se and Sp default to the 'low' setting
                          sd_se_logit = 0.1, sd_sp_logit = 0.1, # SD of Se and Sp default to 'low' setting
                          rho = -0.5,  # Default correlation between Se and Sp
                          floor_prob = 0.5,  # Ensure Se and Sp not below 0.5
                          ceiling_prob = 0.999 # Ensure Se and Sp not perfect
) {

  mu = c(qlogis(mean_se), qlogis(mean_sp))

  # Correlation matrix
  R = matrix(c(1, rho,
                rho, 1), nrow = 2)

  # SDs
  sds = c(sd_se_logit, sd_sp_logit)

  # Covariance matrix
  Sigma = psych::cor2cov(R, sds)

  etas = MASS::mvrnorm(n = di, mu = mu, Sigma = Sigma)
  colnames(etas) = c("eta_se", "eta_sp")

  se_base = plogis(etas[, "eta_se"])
  sp_base = plogis(etas[, "eta_sp"])

  # Keep Se and Sp at or above floor_prob
  se_base = clamp_probability(se_base, floor_prob, ceiling_prob)
  sp_base = clamp_probability(sp_base, floor_prob, ceiling_prob)

  data.frame(
    measure_id = seq_len(di),
    eta_se_base = qlogis(se_base),  # re-logit after bounding so downstream jitter is consistent
    eta_sp_base = qlogis(sp_base),
    se_base = se_base,
    sp_base = sp_base
  )
}


## Study Level Jitter ------------------------------------------------------

jitter_study_accuracy = function(measures, measure_idx,
                                  omega_se, omega_sp,
                                  floor_prob = 0.5, ceiling_prob = 0.999) {
  eta_se = rnorm(length(measure_idx),
                  mean = measures$eta_se_base[measure_idx],
                  sd   = omega_se)
  eta_sp = rnorm(length(measure_idx),
                  mean = measures$eta_sp_base[measure_idx],
                  sd   = omega_sp)

  Se = plogis(eta_se)
  Sp = plogis(eta_sp)

  Se = clamp_probability(Se, floor_prob, ceiling_prob)
  Sp = clamp_probability(Sp, floor_prob, ceiling_prob)

  list(Se = Se, Sp = Sp)
}

## Generate Study Sample Sizes ---------------------------------------------

# Returns sample sizes for the study, using lognormal
sample_sizes_logn = function(sn,  # Number of studies
                              meanlog = 5.50, sdlog = 0.55,  # Default lognormal parameters
                              add_const = 30  # Default constant (i.e., minimum sample size)
) {
  n_raw = rlnorm(sn, meanlog = meanlog, sdlog = sdlog) + add_const
  n = round(n_raw)
  return(n)
}

## Study Level Prevalences -------------------------------------------------

simulate_true_prevalence = function(sn,  # Number of studies
                                     mu_prevalence = 0.015,  # "True" prevalence in a typical study default to 'very low'
                                     tau_logit = 0.5  # Study level heterogeneity in logit space set default to 'low'
) {
  return(plogis(rnorm(sn, mean = qlogis(mu_prevalence), sd = tau_logit)))
}


## Forward Misclassification and Binomial Sampling -------------------------

# Returns list with p_obs and y
forward_misclassification = function(n_i,  # Current sample size
                                      theta_i,  # "True" prevalence
                                      Se_i,  # Current Se
                                      Sp_i   # Current Sp
) {
  p_obs = Se_i * theta_i + (1 - Sp_i) * (1 - theta_i)
  y     = rbinom(length(n_i), size = n_i, prob = p_obs)
  list(p_obs = p_obs, y = y)
}


## Simulate One Data Set ---------------------------------------------------

simulate_dataset = function(
    sn = 20,  # number of studies

    # Sample-size parameters
    meanlog_n = 5.50, sdlog_n = 0.55, add_const_n = 30,

    # Screening parameters
    di = 10,
    mean_se = 0.80, mean_sp = 0.75,
    sd_se_logit = 0.1, sd_sp_logit = 0.1, rho = -0.5,

    # Study-level jitter
    omega_se = 0.2, omega_sp = 0.2,

    # Prevalence
    mu_prevalence = 0.015,            # probability scale
    tau_logit = 0.5,                 # heterogeneity on logit scale

    # Bias in priors and effective sample size of the prior
    sigma_bias = 0, kappa_prior = 200,


    # Bounds for Se/Sp
    floor_prob = 0.5, ceiling_prob = 0.999
) {
  # Generate sample sizes
  n_i = sample_sizes_logn(sn, meanlog_n, sdlog_n, add_const_n)

  # Draw measure-level operating points and assign measures to studies (with replacement)
  measures = draw_measures(di, mean_se, mean_sp, sd_se_logit, sd_sp_logit, rho,
                            floor_prob, ceiling_prob)
  measure_idx = sample(measures$measure_id, size = sn, replace = TRUE)

  acc = jitter_study_accuracy(measures, measure_idx,
                               omega_se = omega_se, omega_sp = omega_sp,
                               floor_prob = floor_prob, ceiling_prob = ceiling_prob)
  Se_i = acc$Se
  Sp_i = acc$Sp

  # True prevalence for each study
  theta_i = simulate_true_prevalence(sn, mu_prevalence, tau_logit)

  # Forward misclassification model
  fm = forward_misclassification(n_i, theta_i, Se_i, Sp_i)

  # data.frame
  data.frame(
    sn = sn,
    di = di,
    mean_se = mean_se,
    mean_sp = mean_sp,
    sd_se_logit = sd_se_logit,
    sd_sp_logit = sd_sp_logit,
    rho = rho,
    omega_se = omega_se,
    omega_sp = omega_sp,
    sigma_bias = sigma_bias,
    kappa_prior = kappa_prior,
    mu_prevalence = mu_prevalence,
    tau_logit = tau_logit,
    study      = seq_len(sn),
    n          = n_i,
    y          = fm$y,
    p_obs      = fm$p_obs,
    theta_true = theta_i,
    Se         = Se_i,
    Sp         = Sp_i,
    measure_id = measure_idx,
    se_base    = measures$se_base[measure_idx],
    sp_base    = measures$sp_base[measure_idx],
    eta_se_base= measures$eta_se_base[measure_idx],
    eta_sp_base= measures$eta_sp_base[measure_idx],
    stringsAsFactors = FALSE
  )
}

## Create Scenarios --------------------------------------------------------

# Map design factors to numeric values
scenario_grid = function(
    sn_levels      = 20,  # Called k in the paper
    accuracy      = 0.80, # Set for Mu Se; Mu Sp is always .05 lower
    sigma_levels  = 0.1,
    omega_levels  = 0.2,
    mu_prev_vals  = 0.015,
    tau = 0.5,
    rho = -0.5,           # The same rho is used for Se/Sp and for bias
    sigma_bias = 0.0,
    kappa_prior = 200
) {
  expand.grid(
    sn         = sn_levels,
    acc_level = accuracy,
    sig_level = sigma_levels,
    omg_level = omega_levels,
    mu_prev   = mu_prev_vals,
    tau = tau,
    rho = rho,
    sigma_bias = sigma_bias,
    kappa_prior = kappa_prior
  )
}


## Create Simulations ------------------------------------------------------

#
# Simulate ALL scenarios (does not fit models)
#
simulate_all_scenarios = function(R = 200, # Number of replications
                                   base_seed = 20250903,  # fixed seed for reproducibility
                                   meanlog_n = 5.50, sdlog_n = 0.55, add_const_n = 30,
                                   di = 10, # Number of measures
                                   floor_prob = 0.5,
                                   ceiling_prob = 0.999,
                                   grid = scenario_grid()) {

  n_scen = nrow(grid)

  out = vector("list", length = n_scen)
  for (s in seq_len(n_scen)) {
    message(sprintf("Scenario %d / %d", s, n_scen))
    sn = grid$sn[s]
    acctag = grid$acc_level[s]
    sigtag = grid$sig_level[s]
    omgtag = grid$omg_level[s]
    mu_p   = grid$mu_prev[s]
    tau_logit = grid$tau[s]
    rho = grid$rho[s]
    sigma_bias = grid$sigma_bias[s]
    kappa_prior = grid$kappa_prior[s]

    mu_se = acctag
    mu_sp = acctag - 0.05 # Note that Sp is always .05 less than Se

    sd_se = sigtag
    sd_sp = sigtag

    om_se = omgtag
    om_sp = omgtag

    # Replicates for this scenario
    datasets = vector("list", length = R)
    for (r in seq_len(R)) {
      set.seed(base_seed + (s - 1) * 1e6 + r)  # reproducible per-scenario, per-replicate
      datasets[[r]] = simulate_dataset(
        sn = sn,
        meanlog_n = meanlog_n, sdlog_n = sdlog_n, add_const_n = add_const_n,
        di = di,
        mean_se = mu_se, mean_sp = mu_sp,
        sd_se_logit = sd_se, sd_sp_logit = sd_sp, rho = rho,
        omega_se = om_se, omega_sp = om_sp,
        mu_prevalence = mu_p,
        tau_logit = tau_logit,
        sigma_bias = sigma_bias, kappa_prior = kappa_prior,
        floor_prob = floor_prob, ceiling_prob = ceiling_prob
      )
    }

    out[[s]] = list(
      scenario_id = s,
      factors = list(
        sn = sn,
        accuracy = acctag,
        sigma_level = sigtag,
        omega_level = omgtag,
        mu_prev = mu_p,
        rho = rho
      ),
      replicates = datasets
    )
  }

  # study        : Study index (1 to sn)
  # n            : Sample size for that study, generated from a log-normal distribution + constant
  # y            : Observed number of "cases" classified positive by the screening tool
  #                y_i ~ Binomial(n_i, p_obs)
  # p_obs        : Expected probability of classification as a case, given misclassification
  #                p_obs = Se_i * theta_true + (1 - Sp_i) * (1 - theta_true)

  # theta_true   : True prevalence in the study (probability scale),
  #                generated from logit(theta_i) = beta0 + u_i,
  #                with u_i ~ Normal(0, tau^2)

  # Se           : Study-specific sensitivity (probability the test correctly detects a true case),
  #                derived from the assigned measure + study-level jitter
  # Sp           : Study-specific specificity (probability the test correctly rejects a non-case),
  #                derived from the assigned measure + study-level jitter

  # measure_id   : Identifier for which screening measure (instrument + cutoff) was assigned to the study
  # se_base      : Baseline sensitivity of the assigned measure before study-level jitter
  # sp_base      : Baseline specificity of the assigned measure before study-level jitter

  # eta_se_base  : Logit of se_base (latent mean on the logit scale, used for generating study-level jitter)
  # eta_sp_base  : Logit of sp_base (latent mean on the logit scale, used for generating study-level jitter)
  out
}


## Generate Univariate Priors for SeSp -------------------------------------

prior_logit_for_study = function(p, # Prior centred on this value
                                 i, # Study ID
                                 nlpar = c("Se", "Sp"), # Is this Se or Sp?
                                 n = 200, # Prior certainty (delta-method approximation to a Beta with this effective sample size)
                                 bounded = FALSE
) {
  nlpar = match.arg(nlpar)
  stopifnot(is.numeric(p), length(p) == 1L, p > 0, p < 1)
  stopifnot(is.numeric(i), length(i) == 1L, i >= 1)

  if (bounded) {
    t        = map_to_inner_t(p)
    center   = qlogis(t)
    sd_logit = sqrt((1 + t) / ((n + 1) * t^2 * (1 - t)))
  } else {
    center   = qlogis(p)
    sd_logit = sqrt(1 / ((n + 1) * p * (1 - p)))
  }

  set_prior(
    sprintf("normal(%0.6f, %0.6f)", center, sd_logit),
    nlpar = nlpar,
    class  = "b",
    coef   = sprintf("study%d", i)
  )
}

prior_logit_for_measure = function(p, # Prior centred on this value
                                   i, # Measure ID
                                   nlpar = c("Se", "Sp"), # Is this Se or Sp?
                                   n = 200, # Prior certainty (delta-method approximation to a Beta with this effective sample size)
                                   bounded = FALSE # TRUE for the 0.5 + 0.5*inv_logit() parameterization
) {
  nlpar = match.arg(nlpar)
  stopifnot(is.numeric(p), length(p) == 1L, p > 0, p < 1)
  stopifnot(is.numeric(i), length(i) == 1L, i >= 1)

  if (bounded) {
    # Se = 0.5 + 0.5*inv_logit(eta)  =>  eta = logit(2p - 1).
    # Delta-method SD of eta for a Beta(mean = p, ess = n) belief on the OUTER p,
    # propagated through the FULL map p -> t -> logit(t):
    t        = map_to_inner_t(p)
    center   = qlogis(t)
    sd_logit = sqrt((1 + t) / ((n + 1) * t^2 * (1 - t)))
  } else {
    # Unbounded: eta = logit(Se); delta-method SD of a Beta(mean = p, ess = n).
    center   = qlogis(p)
    sd_logit = sqrt(1 / ((n + 1) * p * (1 - p)))
  }

  set_prior(
    sprintf("normal(%0.6f, %0.6f)", center, sd_logit),
    nlpar = nlpar,
    class  = "b",
    coef   = sprintf("measure_id%d", i)
  )
}


## Access population prevalence (probability space) and tau from any model ---------------------

# M1-3: plogis(b_Intercept); M4-6,8: plogis(b_pi_Intercept)
.prev_draws = function(fit) {
  dr = as_draws_df(fit)
  if ("b_pi_Intercept" %in% names(dr)) {
    plogis(dr$b_pi_Intercept)
  } else if ("b_Intercept" %in% names(dr)) {
    plogis(dr$b_Intercept)
  } else {
    stop("Could not find population intercept in draws.")
  }
}

# Extract draws for tau
.tau_draws = function(fit) {
  stopifnot(inherits(fit, "brmsfit"))
  dr  = posterior::as_draws_df(fit)
  nms = names(dr)

  # Try the common names in order of specificity
  cand = c(
    "sd_study__pi_Intercept",                 # nlpar = "pi" random intercept by study
    "sd_study__pi__Intercept",                # alternative brms naming
    "sd_study__Intercept",                     # back-up
    "sd_es_id__pi_Intercept", 
    "sd_es_id__pi__Intercept",
    "sd_es_id__Intercept"
  )

  hit = cand[cand %in% nms]
  if (length(hit) == 1L) return(dr[[hit]])

  # Last-resort regex in case of minor naming variants
  hit2 = grep("^sd_study__.*pi.*Intercept$", nms, value = TRUE)
  if (length(hit2) == 1L) return(dr[[hit2]])

  stop(
    "Could not find tau (study-level SD on pi). Checked: ",
    paste(cand, collapse = ", "),
    " and a regex fallback.\n",
    "Available names include: ", paste(utils::head(nms, 40), collapse = ", "),
    "\nTip: inspect names(as_draws_df(fit)) to verify the sd_study__pi_Intercept name."
  )
}


## Access simulation metrics -----------------------------------------------

# Simulation metrics
simulation_metrics = function(fit, dat_raw, y_cor = FALSE) {
  # Sim parameters
  sim_param = dat_raw %>% select(sn:tau_logit) %>% distinct()

  # Prevalence posterior
  p = .prev_draws(fit)

  # Tau posterior
  tau = .tau_draws(fit)

  # Posterior predictive diagnostics
  yrep = posterior_predict(fit, ndraws = 500)   # matrix: draws × observations

  # observed data (vector)
  if(y_cor)
  {
    y_obs = fit$data$cor_y
  } else {
    y_obs = fit$data$y
  }

  tibble(
    # Prevalence metrics
    prev_mean   = mean(p),
    prev_sd     = sd(p),
    prev_lb   = quantile(p, 0.025),
    prev_ub  = quantile(p, 0.975),
    bias = mean(p) - sim_param$mu_prevalence,
    rel_bias = (bias) / sim_param$mu_prevalence,
    sqerror = (bias)**2,
    coverage = prev_lb < sim_param$mu_prevalence & prev_ub > sim_param$mu_prevalence,

    # Heterogeneity metrics
    tau_mean = mean(tau),
    tau_sd = sd(tau),
    tau_lb   = quantile(tau, 0.025),
    tau_ub  = quantile(tau, 0.975),
    tau_bias = mean(tau) - sim_param$tau_logit,
    tau_rel_bias = (tau_bias) / sim_param$tau_logit,
    tau_sqerror = (tau_bias)**2,
    tau_coverage = tau_lb < sim_param$tau_logit & tau_ub > sim_param$tau_logit,

    # Sampler Diagnostics
    maxrhat = max(rhat(fit), na.rm = TRUE),
    proprhat = mean(rhat(fit) > 1.01, na.rm = TRUE),
    minbulk = min(posterior::summarise_draws(posterior::as_draws(fit), "ess_bulk")$ess_bulk, na.rm = TRUE),
    mintail = min(posterior::summarise_draws(posterior::as_draws(fit), "ess_tail")$ess_tail, na.rm = TRUE),
    ndiv = sum(do.call(rbind, rstan::get_sampler_params(fit$fit, inc_warmup = FALSE))[,"divergent__"]),
    maxtreedep = max(do.call(rbind, rstan::get_sampler_params(fit$fit, inc_warmup = FALSE))[,"treedepth__"]),

    # Predictive diagnostics
    ppc_md = mean(rowMeans(yrep)) - mean(y_obs),
    ppc_var_ratio = mean(apply(yrep, 1, var)) / var(y_obs)
  )
}

## Run a single replication of the simulation ------------------------------------------------------

perturb_sesp_study = function(SeSp, sigma_SeSp, eps = 1e-6) {
  # If sigma_SeSp is 0, no perturbation applied, otherwise add random noise
  # based on on a logit transform; latch it on to the eps provided, so
  # Se/Sp doesn't go above 1 or below .5 (where it would be unstable)
  if(sigma_SeSp > 0)
  {
    SeSp_z = qlogis(SeSp) + rnorm(length(SeSp), mean = 0, sd = sigma_SeSp)
    SeSp = plogis(SeSp_z)
    SeSp = clamp_probability(SeSp, .5 + eps, 1 - eps)
  }

  return(SeSp)
}

perturb_sesp_by_measure = function(p, measure_id, sigma, lb = 0.5, eps = 1e-6) {
  # If sigma is 0, don't apply any correction
  if (sigma <= 0) return(p)

  # If sigma is more than 0, generate separate modifers (via a logit link) for each study and
  # then apply them as with the SeSp by study correction, if it goes above or below 1 or .5,
  # replace with that value.
  key = as.character(measure_id)
  lev = unique(key)
  delta = rnorm(length(lev), mean = 0, sd = sigma)  # one draw per measure
  names(delta) = lev

  p_new = plogis(qlogis(p) + delta[key])            # shared shift within measure
  p_new = clamp_probability(p_new, lb + eps, 1 - eps)               # clamp to (lb, 1)
  names(p_new) = NULL
  return(p_new)
}

fit_all_models_one_rep = function(dat_raw, iter = 6000, adapt_delta = 0.999,
                                  is_bounded=FALSE, max_treedepth=20,
                                  stepsize = NULL, stan_m7 = NULL) {
  true_pr = dat_raw$mu_prevalence[1]
  kappa_pr = dat_raw$kappa_prior[1]
  c_sigma_bias = dat_raw$sigma_bias[1]
  common_priors = make_common_priors(qlogis(true_pr))
  dat_raw$is_bounded = is_bounded

  if(is.null(stepsize))
  {
    ctrl_parameters = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  } else {
    ctrl_parameters = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                           step_size = stepsize)
  }

  # Naive corrections
  dat_m1 = dat_raw %>% select(y, n, study)

  # Perturb the per study SeSp if sigma_bias > 0 and correct by study
  dat_m2 = dat_raw %>%
    mutate(Se_p = perturb_sesp_study(Se, c_sigma_bias), Sp_p = perturb_sesp_study(Sp, c_sigma_bias)) %>%
    mutate(cor_y = round(n * rogan_gladen(y/n, Se_p, Sp_p))) %>%
    select(cor_y, n, study) %>% mutate(study = factor(study))

  # Perturb the per measure SeSp if sigma_bias > 0 and correct by measure
  dat_m3 = dat_raw %>%
    mutate(Se_p = perturb_sesp_by_measure(se_base, measure_id, c_sigma_bias), Sp_p = perturb_sesp_by_measure(sp_base, measure_id, c_sigma_bias)) %>%
    mutate(cor_y = round(n * rogan_gladen(y/n, Se_p, Sp_p))) %>%
    select(cor_y, n, study) %>%
    mutate(study = factor(study))

  # Misclassification models
  dat_m4 = dat_raw %>% select(y, n, study, Se, Sp)
  dat_m5 = dat_raw %>% select(y, n, study, measure_id, se_base, sp_base)
  dat_m6 = dat_raw %>% select(y, n, study, measure_id, se_base, sp_base)

  # Priors that depend on the replicate
  pri4 = NULL
  if(is_bounded) {
    pri4 = build_priors_m4_bounded(dat_m4, common_priors, n_eff = kappa_pr, sigma_bias=c_sigma_bias)
  } else {
    pri4 = build_priors_m4(dat_m4, common_priors, n_eff = kappa_pr, sigma_bias=c_sigma_bias)
  }

  pri5 = NULL
  if(is_bounded) {
    pri5 = build_priors_m5_bounded(dat_m5, common_priors, n_eff = kappa_pr, sigma_bias=c_sigma_bias)
  } else {
    pri5 = build_priors_m5(dat_m5, common_priors, n_eff = kappa_pr, sigma_bias=c_sigma_bias)
  }

  pri6 = NULL
  if(is_bounded) {
    pri6 = build_priors_m6_bounded(dat_m6, common_priors, n_eff = kappa_pr, sigma_bias=c_sigma_bias)
  } else {
    pri6 = build_priors_m6(dat_m6, common_priors, n_eff = kappa_pr, sigma_bias=c_sigma_bias)
  }

  # Necessary to avoid these columns informing the model below as data
  dat_m4 = dat_m4 %>% select(-Se, -Sp)

  #
  # Fit Models
  #

  # Model 1: Naive
  fit1 = brm(
    y | trials(n) ~ 1 + (1 | study),
    data = dat_m1,
    family = binomial(),
    backend='cmdstanr',
    prior = common_priors$m123,
    chains = 4,
    cores = 4,
    iter = iter,
    refresh = 0,
    silent = 2,
    control = ctrl_parameters
  )

  # Model 2: Correction based on study
  fit2 = brm(
    cor_y | trials(n) ~ 1 + (1 | study),
    data = dat_m2,
    family = binomial(),
    backend='cmdstanr',
    prior = common_priors$m123,
    chains = 4,
    cores = 4,
    iter = iter,
    refresh = 0,
    silent = 2,
    control = ctrl_parameters
  )

  # Model 3: Correction based on measure mean
  fit3 = brm(
    cor_y | trials(n) ~ 1 + (1 | study),
    data = dat_m3,
    family = binomial(),
    backend='cmdstanr',
    prior = common_priors$m123,
    chains = 4,
    cores = 4,
    iter = iter,
    refresh = 0,
    silent = 2,
    control = ctrl_parameters
  )

  # Model 4: Se/Sp modelled via priors centred on true study values
  c_bf_m4 = NULL
  if(is_bounded) {
    c_bf_m4 = bf_m4_bounded
  } else {
    c_bf_m4 = bf_m4
  }

  fit4 = brm(
    formula = c_bf_m4,
    data = dat_m4 %>% mutate(study = as.character(study)),
    family = binomial(link='identity'),
    backend='cmdstanr',
    prior = pri4,
    chains = 4,
    cores = 4,
    iter = iter,
    refresh = 0,
    silent = 2,
    control = ctrl_parameters
  )

  # Model 5: Se/Sp modelled via priors centred on true measure values
  c_bf_m5 = NULL
  if(is_bounded) {
    c_bf_m5 = bf_m5_bounded
  } else {
    c_bf_m5 = bf_m5
  }

  fit5 = brm(
    formula = c_bf_m5,
    data = dat_m5 %>% mutate(measure_id = as.character(measure_id)),
    family = binomial(link='identity'),
    backend='cmdstanr',
    prior = pri5,
    chains = 4,
    cores = 4,
    iter = iter,
    refresh = 0,
    silent = 2,
    control = ctrl_parameters
  )

  # Model 6: Se/Sp modelled via priors centred on true study values + random effects for Se/Sp
  c_bf_m6 = NULL
  if(is_bounded) {
    c_bf_m6 = bf_m6_bounded
  } else {
    c_bf_m6 = bf_m6
  }

  fit6 = brm(
    formula = c_bf_m6,
    data = dat_m6 %>% mutate(measure_id = as.character(measure_id)),
    family = binomial(link='identity'),
    backend='cmdstanr',
    prior = pri6,
    chains = 4,
    cores = 4,
    iter = iter,
    refresh = 0,
    silent = 2,
    control = ctrl_parameters
  )

  # Model 7: Same as Model 6, but allowing for Se/Sp priors to be correlated
  fit7_metrics = NULL
  if (!is.null(stan_m7)) {
    stan_dat = build_stan_data_sim(
      dat_raw,
      n_eff      = kappa_pr,
      sigma_bias = c_sigma_bias,
      is_bounded = is_bounded
    )

    stan_args = list(
      data             = stan_dat,
      chains           = 4,
      parallel_chains  = 4,
      iter_sampling    = iter %/% 2,
      iter_warmup      = iter %/% 2,
      adapt_delta      = adapt_delta,
      max_treedepth    = max_treedepth,
      refresh          = 0,
      show_messages    = FALSE,
      show_exceptions  = FALSE
    )
    if (!is.null(stepsize)) {
      stan_args$step_size = stepsize
    }

    fit7 = do.call(stan_m7$sample, stan_args)
    fit7_metrics = simulation_metrics_stan(fit7, dat_raw, stan_dat) %>%
      mutate(model = "M7")
  }

  #
  # Summaries (prevalence)
  #
  k=bind_rows(
    simulation_metrics(fit1, dat_raw)  %>% mutate(model = "M1"),
    simulation_metrics(fit2, dat_raw, y_cor=TRUE)  %>% mutate(model = "M2"),
    simulation_metrics(fit3, dat_raw, y_cor=TRUE)  %>% mutate(model = "M3"),
    simulation_metrics(fit4, dat_raw)  %>% mutate(model = "M4"),
    simulation_metrics(fit5, dat_raw)  %>% mutate(model = "M5"),
    simulation_metrics(fit6, dat_raw)  %>% mutate(model = "M6"),
    fit7_metrics
  ) %>% bind_cols(dat_raw %>% select(sn:tau_logit) %>% distinct())

  return(k)
}

## Main simulation functions -----------------------------------------------------------

run_all_rds = function(sims, dir, iter = 6000, adapt_delta = 0.999,
                        skip_num=0, stan_m7 = NULL, is_bounded = FALSE) {
  for(i in (seq_len(length(sims)-skip_num)+skip_num))
  {
    message(sprintf("Starting i = %d", i))
    c_sim = sims[[i]]

    for(j in seq_len(length(c_sim$replicates)))
    {
      fn = sprintf('%sSim%d_%d.rds', dir, i, j)

      if(file.exists((fn))) {
        message(sprintf("Loading r = %d", j)) # Note it says loading because my old code did load, legacy now, it skips
      } else {
        start_time = Sys.time()
        c_rep = c_sim$replicates[[j]]

        # keep original replicate data as list-column
        summ = fit_all_models_one_rep(c_rep, iter = iter,
                                      adapt_delta = adapt_delta, stan_m7 = stan_m7,
                                      is_bounded = is_bounded)
        temp = summ %>%
          mutate(sim_id = i, rep_id = j) %>%
          select(sim_id, rep_id, model, everything()) %>%
          mutate(data = list(c_rep))  # list-column with the full dataset

        saveRDS(temp, fn)
        end_time = Sys.time()
        message(sprintf("     Completed r = %d in %.02f s", j, as.numeric(end_time - start_time, units = "secs")))
      }
    }
  }
}

#
# Attempt to repair simulations, re-running iterations that ran into trouble
#
repair_cycle = function(sims, iter = 12000, adapt_delta = 0.999,
                        step_size=0.002, dir='simulation_data/', cor_dir='corrected_data/',
                        cycles=3, stan_m7 = NULL, is_bounded = FALSE)
{
  for(i in seq_len(cycles))
  {
    results_dat = read_sim_rds(dir=dir)
    results_dat %>%
      filter(ndiv > 0 | maxrhat > 1.01) %>%
      select(sim_id, rep_id) %>%
      distinct() -> bad_trials

    # Clear old corrected files
    for(j in seq_len(nrow(bad_trials)))
    {
      s = bad_trials$sim_id[j]
      r = bad_trials$rep_id[j]
      bad_fn = sprintf('%sSim%d_%d.rds', cor_dir, s, r)

      if(file.exists(bad_fn))
      {
        file.remove(bad_fn)
      }
    }

    if(nrow(bad_trials) > 0)
    {
      fix_all_rds(bad_trials, sims, iter=iter, adapt_delta=adapt_delta,
                  step_size=step_size, dir=dir, cor_dir=cor_dir,
                  stan_m7 = stan_m7, is_bounded = is_bounded)
    } else {
      return(0)
    }
  }

  results_dat = read_sim_rds(dir=dir)
  results_dat %>%
    filter(ndiv > 0 | maxrhat > 1.01) %>%
    select(sim_id, rep_id) %>%
    distinct() -> bad_trials

  return(nrow(bad_trials))
}

fix_all_rds = function(bad_trials, sims, iter = 12000, adapt_delta = 0.999,
                       step_size=0.002, dir='simulation_data/', cor_dir='corrected_data/',
                       stan_m7 = NULL, is_bounded = FALSE) {
  # Keep the previous result before replacing it with a repaired fit.
  if (!dir.exists(cor_dir) && !dir.create(cor_dir, recursive = TRUE)) {
    stop("Could not create backup directory: ", cor_dir)
  }

  nbad_trials = nrow(bad_trials)
  for(i in seq_len(nbad_trials))
  {
    s = bad_trials$sim_id[i]
    r = bad_trials$rep_id[i]
    message(sprintf("Fixing s, r = %d, %d; i = %d of %d", s, r, i, nbad_trials))
    dat_raw = sims[[s]]$replicates[[r]]

    # Original file
    fn = sprintf('%sSim%d_%d.rds', dir, s, r)
    bad_fn = sprintf('%sSim%d_%d.rds', cor_dir, s, r)

    if(file.exists((bad_fn))) {
      message("    File Exists. Skipping")
      next
    }

    if (!file.rename(fn, bad_fn)) {
      stop("Could not back up ", fn, "; the original result has not been replaced.")
    }

    start_time = Sys.time()

    # keep original replicate data as list-column
    summ = fit_all_models_one_rep(dat_raw, iter = iter,
                                  adapt_delta = adapt_delta, stepsize=step_size,
                                  stan_m7 = stan_m7,
                                  is_bounded = is_bounded)
    temp = summ %>%
      mutate(sim_id = s, rep_id = r) %>%
      select(sim_id, rep_id, model, everything()) %>%
      mutate(data = list(dat_raw))  # list-column with the full dataset

    saveRDS(temp, fn)
    end_time = Sys.time()
    message(sprintf("     Completed s, r = %d, %d in %.2f s", s, r, as.numeric(end_time - start_time, units = "secs")))
  }
}

#
# Reads and combine the simulation rds files
#
read_sim_rds = function(dir) {
  pattern = "^Sim([0-9]{1,3})_([0-9]{1,3})\\.rds$"

  # List candidate files
  files = list.files(dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    warning("No matching .rds files found in directory: ", dir)
    return(tibble::tibble())  # return empty tibble
  }

  # Read and combine
  dfs = lapply(files, readRDS)

  # Safety check: enforce data.frame / tibble
  dfs = lapply(dfs, function(x) {
    if (!is.data.frame(x)) stop("One of the RDS files is not a data.frame")
    x
  })

  out = dplyr::bind_rows(dfs, .id = "source_file")
  out$source_file = basename(files)[as.integer(out$source_file)]

  out
}


# Simulation Processing and Plotting Functions ----------------------------


## Main Processing Function ------------------------------------------------

process_simulation = function(
    dir,
    group_vars,
    model_filter   = c("M1", "M3", "M5", "M6"), # Defaults to the basic feasible models
    summary_type   = "prevalence",
    exclude_models = NULL,
    exclude_bad    = FALSE
) {

  # Read & flag bad replications
  results = dplyr::bind_rows(Filter(NROW, lapply(dir, function(d) {
    tryCatch(read_sim_rds(dir = d), warning = function(w) NULL, error = function(e) NULL)
  })))

  bad_reps = results %>%
    dplyr::filter(ndiv > 0 | maxrhat > 1.01) %>%
    dplyr::select(sim_id, rep_id) %>%
    dplyr::distinct()

  if (exclude_bad) {
    results = results %>%
      dplyr::anti_join(bad_reps, by = c("sim_id", "rep_id"))
  }

  # Filter models
  if (!is.null(exclude_models)) {
    results = results %>%
      dplyr::filter(!model %in% exclude_models)
  }

  results = results %>%
    dplyr::filter(model %in% model_filter)

  # Group & summarize
  results_grouped = results %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, "model"))))

  if (summary_type == "prevalence") {

    summary = results_grouped %>%
      dplyr::summarize(
        rmse     = sqrt(mean(sqerror, na.rm = TRUE)),
        ci_width = mean(prev_ub - prev_lb, na.rm = TRUE),
        coverage = mean(coverage, na.rm = TRUE),
        bias     = mean(bias, na.rm = TRUE),
        tau_cov  = mean(tau_coverage, na.rm = TRUE),
        tau_bias = mean(tau_bias, na.rm = TRUE),
        maxrhat  = mean(maxrhat > 1.01, na.rm = TRUE),
        proprhat = mean(proprhat > 0, na.rm = TRUE),
        div      = mean(ndiv > 0, na.rm = TRUE),
        ppc_md   = mean(ppc_md, na.rm = TRUE),
        n        = dplyr::n(),
        .groups  = "drop"
      )

  } else if (summary_type == "moderator") {

    summary = results_grouped %>%
      dplyr::summarize(
        slope_rmse     = sqrt(mean(slope_sqerror, na.rm = TRUE)),
        slope_bias     = mean(slope_bias, na.rm = TRUE),
        slope_coverage = mean(slope_coverage, na.rm = TRUE),
        slope_power    = mean(slope_power, na.rm = TRUE),
        slope_width    = mean(slope_width, na.rm = TRUE),
        prev_rmse      = sqrt(mean(prev_sqerror, na.rm = TRUE)),
        prev_bias      = mean(prev_bias, na.rm = TRUE),
        prev_coverage  = mean(prev_coverage, na.rm = TRUE),
        maxrhat  = mean(maxrhat > 1.01, na.rm = TRUE),
        div      = mean(ndiv > 0, na.rm = TRUE),
        n        = dplyr::n(),
        .groups  = "drop"
      )

  } else {
    stop("summary_type must be 'prevalence' or 'moderator'")
  }

  list(
    summary  = summary,
    bad_reps = bad_reps,
    n_bad    = nrow(bad_reps)
  )
}

make_summary_table <- function(dir,
                               model_filter = "M6",
                               group_vars,
                               average_over = NULL,
                               exclude_bad = TRUE) {

  results <- dplyr::bind_rows(Filter(NROW, lapply(dir, function(d) {
    tryCatch(read_sim_rds(dir = d), warning = function(w) NULL, error = function(e) NULL)
  })))

  if (exclude_bad) {
    bad_reps <- results %>%
      dplyr::filter(ndiv > 0 | maxrhat > 1.01) %>%
      dplyr::distinct(sim_id, rep_id)

    results <- results %>%
      dplyr::anti_join(bad_reps, by = c("sim_id", "rep_id"))
  }

  results <- results %>%
    dplyr::filter(model %in% model_filter)

  results %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      bias     = mean(bias, na.rm = TRUE),
      coverage = mean(coverage, na.rm = TRUE),
      ci_width = mean(prev_ub - prev_lb, na.rm = TRUE),
      n        = dplyr::n(),
      .groups  = "drop"
    )
}

make_moderator_table <- function(dir,
                                 model_filter = c("M1", "Gold", "AllGold", "M6", "M8"),
                                 group_vars   = c("slope_true", "prop_gold", "prior_sp"),
                                 exclude_bad  = TRUE) {

  results <- dplyr::bind_rows(Filter(NROW, lapply(dir, function(d) {
    tryCatch(read_sim_rds(dir = d), warning = function(w) NULL, error = function(e) NULL)
  })))

  if (exclude_bad) {
    bad_reps <- results %>%
      dplyr::filter(ndiv > 0 | maxrhat > 1.01) %>%
      dplyr::distinct(sim_id, rep_id)

    results <- results %>%
      dplyr::anti_join(bad_reps, by = c("sim_id", "rep_id"))
  }

  results <- results %>%
    dplyr::filter(model %in% model_filter)

  results %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, "model")))) %>%
    dplyr::summarise(
      power    = mean(slope_power, na.rm = TRUE),
      coverage = mean(slope_coverage, na.rm = TRUE),
      bias     = mean(slope_bias, na.rm = TRUE),
      width    = mean(slope_width, na.rm = TRUE),
      n        = dplyr::n(),
      .groups  = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from  = model,
      values_from = c(power, coverage, bias, width, n),
      names_glue  = "{model}_{.value}"
    )
}

make_sim8_table <- function(dir = "simulation_data/sim8/",
                            group_vars = c("mean_sp", "mu_prevalence"),
                            exclude_bad = TRUE) {

  results <- read_sim8_rds(dir = dir)

  if (is.null(results) || nrow(results) == 0) {
    stop("No Simulation 8 results found in ", dir)
  }

  if (exclude_bad) {
    results <- results %>%
      dplyr::filter(
        comp_maxrhat  < 1.01, joint_maxrhat < 1.01,
        m6d_maxrhat   < 1.01, m6j_maxrhat   < 1.01,
        comp_ndiv     < 1,    joint_ndiv    < 1,
        m6d_ndiv      < 1,    m6j_ndiv      < 1
      )
  }

  results %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      # M8 prevalence
      m8_prev_bias     = mean(joint_prev_bias, na.rm = TRUE),
      m8_prev_coverage = mean(joint_prev_coverage, na.rm = TRUE),

      # M6j prevalence (M6 with M8-derived priors)
      m6j_prev_bias     = mean(m6j_prev_bias, na.rm = TRUE),
      m6j_prev_coverage = mean(m6j_prev_coverage, na.rm = TRUE),

      # M6d prevalence (M6 with KWGA-derived priors)
      m6d_prev_bias     = mean(m6d_prev_bias, na.rm = TRUE),
      m6d_prev_coverage = mean(m6d_prev_coverage, na.rm = TRUE),

      # KWGA prevalence
      kwga_prev_bias     = mean(kwga_prev_bias, na.rm = TRUE),
      kwga_prev_coverage = mean(kwga_prev_coverage, na.rm = TRUE),

      # Sp recovery
      joint_sp_mae      = mean(abs(joint_sp_bias), na.rm = TRUE),
      joint_sp_bias     = mean(joint_sp_bias, na.rm = TRUE),
      joint_sp_coverage = mean(joint_sp_coverage, na.rm = TRUE),
      disc_sp_bias      = mean(disc_sp_bias_wt, na.rm = TRUE),

      # Se recovery
      joint_se_bias     = mean(joint_se_bias, na.rm = TRUE),
      joint_se_coverage = mean(joint_se_coverage, na.rm = TRUE),

      n = dplyr::n(),
      .groups = "drop"
    )
}

# Prior Generations -------------------------------------------------------

build_priors_m4 = function(dat_fac, priors_common, n_eff = 200, sigma_bias = 0, floor_prob = .5, ceiling_prob = .999, rho_bias = -.5) {
  pri = priors_common$pi

  centers = dat_fac %>%
    dplyr::select(study, Se, Sp) %>%
    dplyr::distinct() %>%
    dplyr::arrange(study)

  n_study = nrow(centers)

  # Bias draws; if sigma is 0, create a 0 matrix instead. Assumes everything is well formed.
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_study, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_study), rep(0, n_study))
  }

  # Compute biased prior centers per study
  eta_se  = stats::qlogis(centers$Se) + E[, 1]
  eta_sp  = stats::qlogis(centers$Sp) + E[, 2]
  se_b    = stats::plogis(eta_se)
  sp_b    = stats::plogis(eta_sp)

  # Optional clamping to match any model constraints (e.g., Se/Sp in [0.5, 1])
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

  for (i in seq_len(n_study)) {
    j = centers$study[i]
    pri = c(
      pri,
      prior_logit_for_study(se_b[i], j, nlpar = "Se", n = n_eff),
      prior_logit_for_study(sp_b[i], j, nlpar = "Sp", n = n_eff)
    )
  }

  pri
}

build_priors_m5 = function(dat_fac,
                            priors_common,
                            n_eff        = 200,   # kappa: effective prior sample size
                            sigma_bias   = 0,     # SD of prior-center bias (logit scale)
                            rho_bias     = -0.5,  # correlation of Se/Sp bias
                            floor_prob   = 0.5,   # clamp external Se/Sp to [floor_prob, ceiling_prob]
                            ceiling_prob = 0.999) {

  pri = priors_common$pi

  centers = dat_fac %>%
    dplyr::select(measure_id, se_base, sp_base) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  n_meas = nrow(centers)

  # Bias draws; if sigma is 0, create a 0 matrix instead. Assumes everything is well formed.
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_meas, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_meas), rep(0, n_meas))
  }

  # Apply bias on logit scale, then back to probability
  eta_se = stats::qlogis(centers$se_base) + E[, 1]
  eta_sp = stats::qlogis(centers$sp_base) + E[, 2]
  se_b   = stats::plogis(eta_se)
  sp_b   = stats::plogis(eta_sp)

  # Clamp to match any model constraints (if used)
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

  # Priors centered at (possibly biased) external Se/Sp for each measure
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

build_priors_m6 = function(dat_fac,
                            priors_common,
                            n_eff        = 200,   # kappa: effective prior sample size
                            sigma_bias   = 0,     # SD of prior-center bias (logit scale)
                            rho_bias     = -0.5,  # correlation of Se/Sp bias
                            floor_prob   = 0.5,   # clamp external Se/Sp to [floor_prob, ceiling_prob]
                            ceiling_prob = 0.999) {

  # baseline priors: prevalence + SD priors for study RE on Se/Sp
  pri = c(priors_common$pi, priors_common$se_sp_sd)

  centers = dat_fac %>%
    dplyr::select(measure_id, se_base, sp_base) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  n_meas = nrow(centers)

  # Bias draws per measure; if sigma is 0, create a 0 matrix instead.
  # Assumes everything is well formed.
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_meas, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_meas), rep(0, n_meas))
  }

  # apply bias on logit scale, then back to probability
  eta_se = stats::qlogis(centers$se_base) + E[, 1]
  eta_sp = stats::qlogis(centers$sp_base) + E[, 2]
  se_b   = stats::plogis(eta_se)
  sp_b   = stats::plogis(eta_sp)

  # optional clamping to match model constraints
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

  # priors on measure-level centers (study-level RE handled by priors_common$se_sp_sd)
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

# Prior Generations for Bounded Models -------------------------------------------------------

# Map external probabilities (>= 0.5) to inner probabilities t ∈ (0, 1):
# t = (p - 0.5) / 0.5 = 2p - 1  (assumes the likelihood uses 0.5 + 0.5*inv_logit(.))
map_to_inner_t = function(p, eps=.0001) {
  clamp_probability(2 * p - 1, eps, 1 - eps)
}

build_priors_m4_bounded = function(dat_fac,
                                    priors_common,
                                    n_eff        = 200,   # kappa: prior effective sample size
                                    sigma_bias   = 0,     # SD of prior-center bias on logit scale
                                    rho_bias     = -0.5,  # correlation of Se/Sp bias
                                    floor_prob   = 0.5,   # clamp external Se/Sp to [floor_prob, ceiling_prob]
                                    ceiling_prob = 0.999) { # keep < 1 to avoid infinite logits


  pri = priors_common$pi

  centers = dat_fac %>%
    dplyr::select(study, Se, Sp) %>%
    dplyr::distinct() %>%
    dplyr::arrange(study)

  n_study = nrow(centers)

  # Draw correlated biases in logit space for Se/Sp; note I'm assuming here
  # that everything is well-formed; if it wasn't, there would be an error, which
  # has never happened; if sigma is 0, create a zeroed out matrix
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_study, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_study), rep(0, n_study))
  }

  # Apply bias to the centers in logit space, then back to probability
  eta_se = stats::qlogis(centers$Se) + E[, 1]
  eta_sp = stats::qlogis(centers$Sp) + E[, 2]
  se_b   = stats::plogis(eta_se)
  sp_b   = stats::plogis(eta_sp)

  # Clamp external probabilities to the max and min
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)


  # Build priors on the inner probabilities (because the NL uses inv_logit on these)
  for (i in seq_len(n_study)) {
    j = centers$study[i]
    pri = c(
      pri,
      prior_logit_for_study(se_b[i], j, nlpar = "Se", n = n_eff, bounded=TRUE),
      prior_logit_for_study(sp_b[i], j, nlpar = "Sp", n = n_eff, bounded=TRUE)
    )
  }

  pri
}

build_priors_m5_bounded = function(dat_fac,
                                    priors_common,
                                    n_eff        = 200,   # kappa: effective prior sample size
                                    sigma_bias   = 0,     # SD of prior-center bias (logit scale)
                                    rho_bias     = -0.5,  # correlation of Se/Sp bias
                                    floor_prob   = 0.5,   # clamp external Se/Sp to [floor_prob, ceiling_prob]
                                    ceiling_prob = 0.999) { # keep < 1 to avoid infinite logits

  pri = priors_common$pi

  centers = dat_fac %>%
    dplyr::select(measure_id, se_base, sp_base) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  n_meas = nrow(centers)

  # Correlated bias draws in logit space (one per measure); this again
  # assumes everything is well-formed. If no pose, create a zeroed
  # out matrix
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_meas, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_meas), rep(0, n_meas))
  }

  # Apply bias on logit scale, then back to probability
  eta_se = stats::qlogis(centers$se_base) + E[, 1]
  eta_sp = stats::qlogis(centers$sp_base) + E[, 2]
  se_b   = stats::plogis(eta_se)
  sp_b   = stats::plogis(eta_sp)

  # Clamp external probabilities
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

  # Priors on the inner probs per measure
  for (i in seq_len(n_meas)) {
    m = centers$measure_id[i]
    pri = c(
      pri,
      prior_logit_for_measure(se_b[i], m, nlpar = "Se", n = n_eff, bounded=TRUE),
      prior_logit_for_measure(sp_b[i], m, nlpar = "Sp", n = n_eff, bounded=TRUE)
    )
  }

  pri
}

build_priors_m6_bounded = function(dat_fac,
                                    priors_common,
                                    n_eff        = 200,   # kappa: effective prior sample size
                                    sigma_bias   = 0,     # SD of prior-center bias (logit scale)
                                    rho_bias     = -0.5,  # correlation of Se/Sp bias
                                    floor_prob   = 0.5,   # clamp external Se/Sp to [floor_prob, ceiling_prob]
                                    ceiling_prob = 0.999) { # keep < 1 to avoid infinite logits


  # baseline priors: prevalence + SD priors for study RE on Se/Sp
  pri = c(priors_common$pi, priors_common$se_sp_sd)

  centers = dat_fac %>%
    dplyr::select(measure_id, se_base, sp_base) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  n_meas = nrow(centers)

  # correlated bias draws in logit space per measure; if sigma is 0
  # create a zero matrix
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_meas, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_meas), rep(0, n_meas))
  }

  # apply bias on logit scale, then back to probability
  eta_se = stats::qlogis(centers$se_base) + E[, 1]
  eta_sp = stats::qlogis(centers$sp_base) + E[, 2]
  se_b   = stats::plogis(eta_se)
  sp_b   = stats::plogis(eta_sp)

  # clamp
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

  # priors on inner probs per measure
  for (i in seq_len(n_meas)) {
    m = centers$measure_id[i]
    pri = c(
      pri,
      prior_logit_for_measure(se_b[i], m, nlpar = "Se", n = n_eff, bounded=TRUE),
      prior_logit_for_measure(sp_b[i], m, nlpar = "Sp", n = n_eff, bounded=TRUE)
    )
  }

  pri
}

## Create common priors used across models ---------------------------------

make_common_priors = function(true_pr) {
  list(
    # M1–M3 (pi only); priors are centred on the true value, see paper for justification
    # this creates priors similar to those used in our past work.
    m123 = c(
      set_prior(sprintf("normal(%f, 1.5)", true_pr), class = "Intercept"),
      prior(normal(0, 1), class = "sd")
    ),
    # M4–M6: pi mean & sd (nlpar)
    pi = c(
      set_prior(sprintf("normal(%f, 1.5)", true_pr), nlpar = "pi"),
      prior(normal(0, 1), nlpar = "pi", class = "sd")
    ),

    # sd priors for Se/Sp study REs where present
    se_sp_sd = c(
      prior(normal(0, .5), nlpar = "Se", class = "sd"),
      prior(normal(0, .5), nlpar = "Sp", class = "sd")
    )
  )
}

# M7 Helper: Build Stan data from simulation replicate ---------------------

build_stan_data_sim = function(dat_raw,
                                n_eff        = 200,
                                sigma_bias   = 0,
                                rho_bias     = -0.5,
                                is_bounded   = FALSE,
                                floor_prob   = 0.5,
                                ceiling_prob = 0.999) {

  # Unique measures and their base Se/Sp
  centers = dat_raw %>%
    dplyr::select(measure_id, se_base, sp_base) %>%
    dplyr::distinct() %>%
    dplyr::arrange(measure_id)

  n_meas = nrow(centers)

  # Remap measure_id to contiguous 1:M
  measure_map = setNames(seq_len(n_meas), centers$measure_id)
  measure_stan = as.integer(measure_map[as.character(dat_raw$measure_id)])

  # Correlated bias perturbation (see build_priors_m5/m6 for something similar)
  if (sigma_bias > 0) {
    Sigma = matrix(c(1, rho_bias, rho_bias, 1), 2, 2) * sigma_bias^2
    E = MASS::mvrnorm(n = n_meas, mu = c(0, 0), Sigma = Sigma)
  } else {
    E = cbind(rep(0, n_meas), rep(0, n_meas))
  }

  eta_se = qlogis(centers$se_base) + E[, 1]
  eta_sp = qlogis(centers$sp_base) + E[, 2]
  se_b   = plogis(eta_se)
  sp_b   = plogis(eta_sp)

  # Clamp for safety
  se_b = clamp_probability(se_b, floor_prob, ceiling_prob)
  sp_b = clamp_probability(sp_b, floor_prob, ceiling_prob)

  # Compute prior parameters per measure
  pmu_se = psd_se = pmu_sp = psd_sp = numeric(n_meas)

  for (m in seq_len(n_meas)) {
    if (is_bounded) {
      # Map to inner probability: t = 2p - 1
      t_se = map_to_inner_t(se_b[m])
      t_sp = map_to_inner_t(sp_b[m])
      pmu_se[m] = qlogis(t_se)
      pmu_sp[m] = qlogis(t_sp)
      psd_se[m] = sqrt((1 + t_se) / ((n_eff + 1) * t_se^2 * (1 - t_se)))
      psd_sp[m] = sqrt((1 + t_sp) / ((n_eff + 1) * t_sp^2 * (1 - t_sp)))
    } else {
      pmu_se[m] = qlogis(se_b[m])
      pmu_sp[m] = qlogis(sp_b[m])
      psd_se[m] = sqrt(1 / ((n_eff + 1) * se_b[m] * (1 - se_b[m])))
      psd_sp[m] = sqrt(1 / ((n_eff + 1) * sp_b[m] * (1 - sp_b[m])))
    }
  }

  list(
    N              = nrow(dat_raw),
    y              = dat_raw$y,
    n              = dat_raw$n,
    M              = n_meas,
    measure        = measure_stan,
    prior_mu_se    = pmu_se,
    prior_mu_sp    = pmu_sp,
    prior_sd_se    = psd_se,
    prior_sd_sp    = psd_sp,
    prior_mu_pi    = qlogis(dat_raw$mu_prevalence[1]),
    prior_sd_pi    = 1.5,
    prior_mu_z_rho = atanh(-0.5),
    prior_sd_z_rho = 0.3
  )
}

# M7 Helper: Compute simulation metrics from cmdstanr fit ------------------

simulation_metrics_stan = function(fit, dat_raw, stan_data) {
  sim_param = dat_raw %>% dplyr::select(sn:tau_logit) %>% dplyr::distinct()

  # Prevalence posterior: plogis(mu_pi)
  p   = plogis(as.vector(fit$draws(variables = "mu_pi", format = "matrix")))

  # Tau posterior: tau_pi

  tau = as.vector(fit$draws(variables = "tau_pi", format = "matrix"))

  # Posterior predictive via p_apparent draws
  p_app       = fit$draws(variables = "p_apparent", format = "matrix")
  n_draws_ppc = min(nrow(p_app), 500)
  idx         = sample(nrow(p_app), n_draws_ppc)
  yrep        = matrix(NA_real_, n_draws_ppc, ncol(p_app))
  for (d in seq_len(n_draws_ppc)) {
    yrep[d, ] = rbinom(ncol(p_app), size = stan_data$n, prob = p_app[idx[d], ])
  }
  y_obs = stan_data$y

  # Sampler diagnostics
  sd_mat = fit$sampler_diagnostics(format = "matrix")
  summ   = fit$summary()

  tibble::tibble(
    # Prevalence metrics
    prev_mean    = mean(p),
    prev_sd      = sd(p),
    prev_lb      = quantile(p, 0.025),
    prev_ub      = quantile(p, 0.975),
    bias         = mean(p) - sim_param$mu_prevalence,
    rel_bias     = bias / sim_param$mu_prevalence,
    sqerror      = bias^2,
    coverage     = prev_lb < sim_param$mu_prevalence & prev_ub > sim_param$mu_prevalence,

    # Heterogeneity metrics
    tau_mean     = mean(tau),
    tau_sd       = sd(tau),
    tau_lb       = quantile(tau, 0.025),
    tau_ub       = quantile(tau, 0.975),
    tau_bias     = mean(tau) - sim_param$tau_logit,
    tau_rel_bias = tau_bias / sim_param$tau_logit,
    tau_sqerror  = tau_bias^2,
    tau_coverage = tau_lb < sim_param$tau_logit & tau_ub > sim_param$tau_logit,

    # Sampler diagnostics
    maxrhat      = max(summ$rhat, na.rm = TRUE),
    proprhat     = mean(summ$rhat > 1.01, na.rm = TRUE),
    minbulk      = min(summ$ess_bulk, na.rm = TRUE),
    mintail      = min(summ$ess_tail, na.rm = TRUE),
    ndiv         = sum(sd_mat[, "divergent__"]),
    maxtreedep   = max(sd_mat[, "treedepth__"]),

    # Predictive diagnostics
    ppc_md        = mean(rowMeans(yrep)) - mean(y_obs),
    ppc_var_ratio = mean(apply(yrep, 1, var)) / var(y_obs)
  )
}

# ---- Helper: find gold/screener columns in draws --------------------------

# Detect column names for gold and screener intercepts in posterior draws.
# Works for both nonlinear models (b_pi_is_goldTRUE) and standard
# binomial models (b_is_goldTRUE).
find_gold_screen_cols <- function(nms) {

  gold_candidates <- c(
    # Nonlinear (nl) models: pi submodel
    "b_pi_is_goldTRUE", "b_pi_is_gold1", "b_pi_is_goldtrue",
    "b_pi_is_gold.TRUE.", "b_pi_is_goldTRUE.",
    # Standard binomial models
    "b_is_goldTRUE", "b_is_gold1", "b_is_goldtrue",
    "b_is_gold.TRUE.", "b_is_goldTRUE."
  )

  screen_candidates <- c(
    # Nonlinear (nl) models
    "b_pi_is_goldFALSE", "b_pi_is_gold0", "b_pi_is_goldfalse",
    "b_pi_is_gold.FALSE.", "b_pi_is_goldFALSE.",
    # Standard binomial models
    "b_is_goldFALSE", "b_is_gold0", "b_is_goldfalse",
    "b_is_gold.FALSE.", "b_is_goldFALSE."
  )

  gold_col <- intersect(gold_candidates, nms)[1]
  screen_col <- intersect(screen_candidates, nms)[1]

  # Fallback: grep for anything containing is_gold
  if (is.na(gold_col)) {
    gold_col <- grep("is_gold.*TRUE|is_gold.*1$", nms, value = TRUE, ignore.case = TRUE)[1]
  }
  if (is.na(screen_col)) {
    screen_col <- grep("is_gold.*FALSE|is_gold.*0$", nms, value = TRUE, ignore.case = TRUE)[1]
  }

  if (is.na(gold_col) || is.na(screen_col)) {
    available <- grep("is_gold|gold|screen", nms, value = TRUE, ignore.case = TRUE)
    stop("Could not find gold/screener intercepts.\n  Available columns matching 'gold'/'screen': ",
         if (length(available) > 0) paste(available, collapse = ", ") else "(none)",
         "\n  All 'b_' columns: ",
         paste(grep("^b_", nms, value = TRUE), collapse = ", "))
  }

  list(gold = gold_col, screen = screen_col)
}


# Helper: Clamp for Proportions and Safe qlogis ---------------------------

clamp_probability = function(p, lower = NULL, upper = NULL) {
  if (!is.null(lower))
    p = pmax(p, lower)
  if (!is.null(upper))
    p = pmin(p, upper)
  p
}

safe_qlogis = function(p, eps = 1e-4) {
  qlogis(clamp_probability(p, eps, 1 - eps))
}


# Models ------------------------------------------------------------------

bf_m4 = bf(
  y | trials(n) ~ (inv_logit(pi) * inv_logit(Se)) + ((1 - inv_logit(pi)) * (1 - inv_logit(Sp))),
  pi ~ 1 + (1|study),
  Se ~ study - 1,
  Sp ~ study - 1,
  nl = TRUE
)

bf_m5 = bf(
  y | trials(n) ~ (inv_logit(pi) * inv_logit(Se)) + ((1 - inv_logit(pi)) * (1 - inv_logit(Sp))),
  pi ~ 1 + (1|study),
  Se ~ measure_id - 1,
  Sp ~ measure_id - 1,
  nl = TRUE
)

bf_m6 = bf(
  y | trials(n) ~ (inv_logit(pi) * inv_logit(Se)) + ((1 - inv_logit(pi)) * (1 - inv_logit(Sp))),
  pi ~ 1 + (1 | s | study),
  Se ~ measure_id - 1 + (1 | s | study),
  Sp ~ measure_id - 1 + (1 | s | study),
  nl = TRUE
)

# Not used, this is for a model that assumes known Se and Sp. Here just in case
# someone wishes to know how to do this:
#
# bf_known = bf(
#   y | trials(n) ~ (inv_logit(pi) * Se) + ((1 - inv_logit(pi)) * (1 - Sp)),
#   pi ~ 1 + (1|study),
#   nl = TRUE
# )
#
# fit_known = brm(
#   formula = bf_known,
#   data = dat_known,
#   family = binomial(link='identity'),
#   backend='cmdstanr',
#   prior = priors_common$pi, # See how these are generated in the simulation code
#   chains = 4,
#   cores = 4,
#   iter = iter,  # See simulation code
#   refresh = 0,
#   silent = 2,
#   control = ctrl_parameters # See simulation code
# )

# Models w/ Support at 0.5 ------------------------------------------------------------------

bf_m4_bounded = bf(
  y | trials(n) ~ (inv_logit(pi) * (0.5 + 0.5 * inv_logit(Se))) +
    ((1 - inv_logit(pi)) * (1 - (0.5 + 0.5 * inv_logit(Sp)))),
  pi ~ 1 + (1 | study),
  Se ~ study - 1,
  Sp ~ study - 1,
  nl = TRUE
)

bf_m5_bounded = bf(
  y | trials(n) ~ (inv_logit(pi) * (0.5 + 0.5 * inv_logit(Se))) +
    ((1 - inv_logit(pi)) * (1 - (0.5 + 0.5 * inv_logit(Sp)))),
  pi ~ 1 + (1|study),
  Se ~ measure_id - 1,
  Sp ~ measure_id - 1,
  nl = TRUE
)

bf_m6_bounded = bf(
  y | trials(n) ~ (inv_logit(pi) * (0.5 + 0.5 * inv_logit(Se))) +
    ((1 - inv_logit(pi)) * (1 - (0.5 + 0.5 * inv_logit(Sp)))),
  pi ~ 1 + (1 | s | study),
  Se ~ measure_id - 1 + (1 | s | study),
  Sp ~ measure_id - 1 + (1 | s | study),
  nl = TRUE
)

# Stan Models -------------------------------------------------------------

# These are extra models not reported in the paper, which use multivariate
# priors on Se and Sp. The model performed no better than M6 where tested.
stan_m7_bounded   = cmdstan_model("stan/rogan_gladen_bvn_measures_clean.stan")
stan_m7_unbounded = cmdstan_model("stan/rogan_gladen_bvn_measures_unbounded.stan")

# Load Specialized Functions for later simulations and diagnostics ------------------------

# The current file was getting too big, so I broke it into parts
# specific to different simulations
source("scripts/0. Setup/0.1 Functions for Simulation 4 and 5.R")
source("scripts/0. Setup/0.2 Functions for Simulation 6.R")
source("scripts/0. Setup/0.3 Functions for Simulation 7.R")
source("scripts/0. Setup/0.4 Functions for Simulation 8.R")
source("scripts/0. Setup/0.5 Functions for SeSp Diagnostics.R")
