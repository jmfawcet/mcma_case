# ---- KWGA: Weighted average prevalence across Se/Sp pairs -----------------

# Compute the KWGA prevalence estimate from Se/Sp weights
# Two modes (first is just faster):
# 1. "analytic": Use the RG-corrected draws from the comparison model
#      (fast, but doesn't account for M6 hierarchical structure)
# 2. "model": Use prevalence draws from pre-fit corrected models at
#      each Se/Sp pair (slower but propagates uncertainty better)
#
# resample (analytic mode only):
#   "independent" - original implementation I used: sample an Se/Sp pair by weight and,
#                   independently, a posterior draw of the screening prevalence
#   "joint"       - sample (posterior draw, Se/Sp pair) together with probability
#                   proportional to prior weight x kernel score, so each draw is
#                   conditioned on its agreement with the gold studies
#   "auto"        - "joint" when the weights were built with clamp_scoring = FALSE
#                   (the default, see compute_sesp_weights), otherwise "independent"
# weights_df may be the grid (as fit_sim8 passes it) or the full list returned by
# compute_sesp_weights(); bandwidth defaults to the one stored with the weights.
kwga_prevalence = function(weights_df,
                            fit_comp = NULL,
                            corrected_fits = NULL,
                            mode = c("analytic", "model"),
                            n_draws = 4000,
                            resample = c("auto", "independent", "joint"),
                            bandwidth = NULL) {

  mode = match.arg(mode)
  resample = match.arg(resample)

  # Accept either the grid or the full compute_sesp_weights() list
  weights_obj = NULL
  if (!is.data.frame(weights_df) && is.list(weights_df) && !is.null(weights_df$grid)) {
    weights_obj = weights_df
    weights_df = weights_obj$grid
  }
  clamp_scoring = if (!is.null(weights_obj)) weights_obj$clamp_scoring else attr(weights_df, "clamp_scoring")
  if (is.null(clamp_scoring)) clamp_scoring = TRUE
  if (is.null(bandwidth)) {
    bandwidth = if (!is.null(weights_obj)) weights_obj$bandwidth else attr(weights_df, "bandwidth")
  }
  if (resample == "auto") resample = if (isFALSE(clamp_scoring)) "joint" else "independent"

  # Filter to Se/Sp pairs with non-negligible weight
  active = weights_df %>% filter(weight > 1e-6)

  if (nrow(active) == 0) {
    stop("No Se/Sp pairs with weight > 1e-6. Check bandwidth or input data.")
  }

  if (mode == "analytic") {

    # Extract screener (and gold) draws
    p_gold = NULL
    if (inherits(fit_comp, "brmsfit")) {
      draws = as_draws_df(fit_comp)
      nms = names(draws)
      cols = find_gold_screen_cols(nms)
      p_screen = plogis(draws[[cols$screen]])
      p_gold   = plogis(draws[[cols$gold]])
    } else if (is.list(fit_comp)) {
      p_screen = fit_comp$screen_draws
      p_gold   = fit_comp$gold_draws
    } else if (is.null(fit_comp) && !is.null(weights_obj)) {
      p_screen = weights_obj$screen_draws
      p_gold   = weights_obj$gold_draws
    } else {
      stop("fit_comp must be a brmsfit or list (or pass the compute_sesp_weights() list as weights_df)")
    }

    if (resample == "independent") {

      # Mixture: sample Se/Sp pairs according to weights, then correct
      pair_idx = sample(seq_len(nrow(active)),
                         size = n_draws,
                         replace = TRUE,
                         prob = active$weight)

      # For each sampled pair, pick a random posterior draw and correct it
      draw_idx = sample(length(p_screen), size = n_draws, replace = TRUE)

      kwga_draws = vapply(seq_len(n_draws), function(k) {
        se = active$se[pair_idx[k]]
        sp = active$sp[pair_idx[k]]
        corrected = (p_screen[draw_idx[k]] + sp - 1) / (se + sp - 1)
        max(0, min(1, corrected))
      }, numeric(1))

    } else {

      # Joint resampling: sample (posterior draw, Se/Sp pair) together with
      # probability proportional to prior weight x kernel score, scored the same
      # way the weights were (clamped or unclamped)
      if (is.null(p_gold) || length(p_gold) != length(p_screen)) {
        stop("Joint resampling needs gold draws paired with the screening draws (same length)")
      }
      if (is.null(bandwidth)) {
        bandwidth = sd(p_gold)
        if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth = 0.02
      }
      n_post = length(p_screen)
      n_act  = nrow(active)
      log_prior = if (is.null(active$log_prior)) rep(0, n_act) else active$log_prior
      raw = outer(p_screen, active$sp - 1, "+") /
        matrix(active$se + active$sp - 1, nrow = n_post, ncol = n_act, byrow = TRUE)
      scored = if (isFALSE(clamp_scoring)) raw else clamp_probability(raw, 0, 1)
      score = dnorm(p_gold - scored, mean = 0, sd = bandwidth, log = TRUE) +
        matrix(log_prior, nrow = n_post, ncol = n_act, byrow = TRUE)
      idx = sample.int(n_post * n_act, size = n_draws, replace = TRUE,
                       prob = exp(score - max(score)))
      kwga_draws = pmin(1, pmax(0, raw[idx]))
    }

  } else {

    stopifnot(!is.null(corrected_fits))

    # Sample Se/Sp pairs according to weights
    pair_idx = sample(seq_len(nrow(active)),
                       size = n_draws,
                       replace = TRUE,
                       prob = active$weight)

    # For each sampled pair, draw from the corresponding corrected model
    kwga_draws = vapply(seq_len(n_draws), function(k) {
      se = active$se[pair_idx[k]]
      sp = active$sp[pair_idx[k]]

      # Find corresponding model (match Se/Sp grid labels)
      label = sprintf("Se%.2f_Sp%.2f", se, sp)
      fit_name = grep(label, names(corrected_fits), value = TRUE, fixed = TRUE)

      if (length(fit_name) == 0) {
        warning(sprintf("No model found for %s", label), call. = FALSE)
        return(NA_real_)
      }

      fit = corrected_fits[[fit_name[1]]]
      dr = as_draws_df(fit)
      prev = plogis(dr$b_pi_Intercept)
      sample(prev, 1)
    }, numeric(1))
  }

  kwga_draws = kwga_draws[!is.na(kwga_draws)]

  if (length(kwga_draws) == 0) {
    stop("All KWGA draws are NA. Check corrected_fits names match grid labels.")
  }

  list(
    summary = tibble(
      mean = mean(kwga_draws),
      median = median(kwga_draws),
      ci_lb = unname(quantile(kwga_draws, 0.025)),
      ci_ub = unname(quantile(kwga_draws, 0.975)),
      sd = sd(kwga_draws),
      n_effective_models = 1 / sum(active$weight^2)
    ),
    draws = kwga_draws,
    weights = weights_df,
    resample = if (mode == "analytic") resample else NA_character_
  )
}

# SeSp Weight Functions ---------------------------------------------------

weighted_quantile = function(x, w, probs = 0.5) {
  o = order(x)
  x = x[o]
  w = w[o]
  cw = cumsum(w) / sum(w)
  approx(cw, x, xout = probs, method = "linear", ties = "ordered")$y
}

weighted_kde_mode = function(x, w, n = 512, adjust = 1) {
  if (length(unique(x)) == 1L) return(unique(x))
  dd = density(
    x,
    weights = w / sum(w),
    n = n,
    adjust = adjust,
    from = min(x),
    to = max(x)
  )
  dd$x[which.max(dd$y)]
}

compute_sesp_weights = function(fit_comp,
                                 se_vals,
                                 sp_vals,
                                 bandwidth = NULL,
                                 prior_weights = NULL,
                                 top_mass = 0.90,
                                 clamp_scoring = FALSE) {

  # Extract posterior draws
  if (inherits(fit_comp, "brmsfit")) {
    draws = as_draws_df(fit_comp)
    nms = names(draws)
    cols = find_gold_screen_cols(nms)

    gold_draws   = plogis(draws[[cols$gold]])
    screen_draws = plogis(draws[[cols$screen]])

  } else if (is.list(fit_comp)) {
    gold_draws   = fit_comp$gold_draws
    screen_draws = fit_comp$screen_draws
  } else {
    stop("fit_comp must be a brmsfit or a list with $gold_draws and $screen_draws")
  }

  gold_draws   = as.numeric(gold_draws)
  screen_draws = as.numeric(screen_draws)

  n_draws = length(gold_draws)
  if (length(screen_draws) != n_draws) {
    stop("gold_draws and screen_draws must have same length")
  }
  if (n_draws < 1) stop("No posterior draws found")

  # Set bandwidth to be used with the normal kernel down below; this is meant to help tune
  # how much mismatch is tolerated between the corrections and the gold standard values.
  # If no bandwidth is given, we use the sd of the gold draws, which is to say less certainty in
  # the gold standard will result in greater tolerance. This seemed reasonable to me, but
  # it's possible there is a better default (or want to encourage reasonable behaviour)
  if (is.null(bandwidth)) {
    bandwidth = sd(gold_draws)
    if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth = 0.02
  }

  # Build grid
  grid = expand.grid(
    se = se_vals,
    sp = sp_vals,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    filter(se + sp > 1)

  n_grid = nrow(grid)
  if (n_grid == 0) stop("No valid (se, sp) pairs with se + sp > 1.")

  # Prior weights on the gride (not used, but a good option to implement)
  if (is.null(prior_weights)) {
    grid$log_prior = rep(0, n_grid)
  } else if (is.vector(prior_weights) && length(prior_weights) == n_grid) {
    if (any(prior_weights < 0) || sum(prior_weights) <= 0) {
      stop("prior_weights must be nonnegative and sum to > 0")
    }
    pw = prior_weights / sum(prior_weights)
    grid$log_prior = log(pw)
  } else {
    warning(sprintf(
      "prior_weights length (%d) != grid size (%d); using uniform",
      length(prior_weights), n_grid
    ))
    grid$log_prior = rep(0, n_grid)
  }

  # Explicit matrices; no sweep()
  denom = grid$se + grid$sp - 1

  screen_mat = matrix(screen_draws, nrow = n_draws, ncol = n_grid)
  sp_mat     = matrix(grid$sp - 1, nrow = n_draws, ncol = n_grid, byrow = TRUE)
  denom_mat  = matrix(denom, nrow = n_draws, ncol = n_grid, byrow = TRUE)
  gold_mat   = matrix(gold_draws, nrow = n_draws, ncol = n_grid, byrow = FALSE)

  corrected_raw = (screen_mat + sp_mat) / denom_mat
  corrected_raw = matrix(corrected_raw, nrow = n_draws, ncol = n_grid)
  corrected = clamp_probability(corrected_raw, 0, 1)
  corrected = matrix(corrected, nrow = n_draws, ncol = n_grid)

  # clamp_scoring = FALSE (default): the unclamped corrected value is scored, so
  # an Se/Sp pair that implies a negative prevalence is penalised by its full
  # distance from the gold draws. clamp_scoring = TRUE scores the clamped value
  # instead, treating such a pair as if it implied zero prevalence. The per-cell
  # summary columns below always use the clamped values.
  discrepancy = gold_mat - (if (isTRUE(clamp_scoring)) corrected else corrected_raw)
  discrepancy = matrix(discrepancy, nrow = n_draws, ncol = n_grid)

  # This is where the kernel from earlier is used, to penalize draws based on
  # distance from 0 within a normal distribution with a set SD
  log_lik_mat = dnorm(discrepancy, mean = 0, sd = bandwidth, log = TRUE)
  log_lik_mat = matrix(log_lik_mat, nrow = n_draws, ncol = n_grid)

  # Stable log-mean-exp by column without colMeans/apply fragility
  max_ll = vapply(seq_len(n_grid), function(k) max(log_lik_mat[, k]), numeric(1))

  log_ml = vapply(seq_len(n_grid), function(k) {
    m = max_ll[k]
    m + log(mean(exp(log_lik_mat[, k] - m)))
  }, numeric(1))

  corrected_mean = vapply(seq_len(n_grid), function(k) {
    mean(corrected[, k])
  }, numeric(1))

  corrected_median = vapply(seq_len(n_grid), function(k) {
    median(corrected[, k])
  }, numeric(1))

  corrected_ci_lb = vapply(seq_len(n_grid), function(k) {
    unname(quantile(corrected[, k], 0.025))
  }, numeric(1))

  corrected_ci_ub = vapply(seq_len(n_grid), function(k) {
    unname(quantile(corrected[, k], 0.975))
  }, numeric(1))

  grid$log_ml = log_ml
  grid$corrected_mean = corrected_mean
  grid$corrected_median = corrected_median
  grid$corrected_ci_lb = corrected_ci_lb
  grid$corrected_ci_ub = corrected_ci_ub

  # Normalize to weights
  log_post = grid$log_ml + grid$log_prior
  max_lp = max(log_post)
  grid$weight = exp(log_post - max_lp)
  grid$weight = grid$weight / sum(grid$weight)

  grid = grid %>%
    arrange(desc(weight)) %>%
    mutate(cumweight = cumsum(weight))

  # Marginal distributions
  se_marg = grid %>%
    group_by(se) %>%
    summarize(w = sum(weight), .groups = "drop") %>%
    arrange(se)

  sp_marg = grid %>%
    group_by(sp) %>%
    summarize(w = sum(weight), .groups = "drop") %>%
    arrange(sp)

  # Estimators
  joint_map = grid %>% slice_max(weight, n = 1, with_ties = FALSE)

  se_map_marg = se_marg %>% slice_max(w, n = 1, with_ties = FALSE) %>% pull(se)
  sp_map_marg = sp_marg %>% slice_max(w, n = 1, with_ties = FALSE) %>% pull(sp)

  se_mean = sum(grid$se * grid$weight)
  sp_mean = sum(grid$sp * grid$weight)

  se_median = weighted_quantile(se_marg$se, se_marg$w, probs = 0.5)
  sp_median = weighted_quantile(sp_marg$sp, sp_marg$w, probs = 0.5)

  se_mode_kde = weighted_kde_mode(se_marg$se, se_marg$w)
  sp_mode_kde = weighted_kde_mode(sp_marg$sp, sp_marg$w)

  se_mean_round = se_vals[which.min(abs(se_vals - se_mean))]
  sp_mean_round = sp_vals[which.min(abs(sp_vals - sp_mean))]

  top_grid = grid %>% filter(cumweight <= top_mass | row_number() == 1)
  top_w = top_grid$weight / sum(top_grid$weight)

  se_top_mean = sum(top_grid$se * top_w)
  sp_top_mean = sum(top_grid$sp * top_w)

  estimates = tibble(
    method = c(
      "joint_map",
      "marginal_map",
      "weighted_mean",
      "weighted_median",
      "kde_mode",
      "weighted_mean_rounded",
      "top_mass_mean"
    ),
    se_est = c(
      joint_map$se,
      se_map_marg,
      se_mean,
      se_median,
      se_mode_kde,
      se_mean_round,
      se_top_mean
    ),
    sp_est = c(
      joint_map$sp,
      sp_map_marg,
      sp_mean,
      sp_median,
      sp_mode_kde,
      sp_mean_round,
      sp_top_mean
    )
  )

  # Carry the settings on the grid too, so kwga_prevalence() can resolve
  # resample = "auto" when it is handed only weights$grid
  attr(grid, "clamp_scoring") = isTRUE(clamp_scoring)
  attr(grid, "bandwidth") = bandwidth

  list(
    grid = grid,
    estimates = estimates,
    gold_draws = gold_draws,
    screen_draws = screen_draws,
    bandwidth = bandwidth,
    clamp_scoring = isTRUE(clamp_scoring)
  )
}
