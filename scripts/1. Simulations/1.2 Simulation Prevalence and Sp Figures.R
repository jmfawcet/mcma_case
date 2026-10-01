# ============================================================================
# Simulation Plotting Functions (All Simulations)
# ============================================================================
#
# Each function:
#   - Loads data via read_sim_rds() internally
#   - Optionally excludes bad replications (maxrhat > 1.01 or ndiv > 0)
#   - Excludes specified models
#   - Incorporates violin plots
#   - Returns a ggplot object
#
# ============================================================================


# ---- Shared helpers ------------------------------------------------------

# Load and filter simulation results
#
# @param dir           Directory containing .rds files
# @param exclude_bad   If TRUE, drop a replication for all models when any model
#   has Rhat > 1.01 or divergences
# @param exclude_models Character vector of models to drop
# @param model_order   Character vector specifying model factor levels
load_sim <- function(dir,
                     exclude_bad    = TRUE,
                     exclude_models = NULL,
                     model_order    = NULL) {
  dat <- read_sim_rds(dir = dir)

  if (exclude_bad) {
    bad_reps <- dat %>%
      dplyr::filter(ndiv > 0 | maxrhat > 1.01) %>%
      dplyr::distinct(sim_id, rep_id)

    dat <- dat %>%
      dplyr::anti_join(bad_reps, by = c("sim_id", "rep_id"))
  }

  if (!is.null(exclude_models)) {
    dat <- dat %>% dplyr::filter(!model %in% exclude_models)
  }

  if (!is.null(model_order)) {
    dat <- dat %>% dplyr::mutate(model = factor(model, levels = model_order))
  }

  dat
}

# Common ggplot theme for all simulation plots
sim_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) %+replace%
    ggplot2::theme(
      panel.border       = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.7),
      panel.grid.minor   = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "black"),
      strip.text         = ggplot2::element_text(color = "white", face = "bold")
    )
}

# Helper to build ordered delta labels
# Returns a factor with levels ordered from negative to positive
ordered_delta_lab <- function(vals, prefix = "\u03B4Sp") {
  labs <- paste0(prefix, " = ", ifelse(vals > 0, "+", ""), sprintf("%.2f", vals))
  lvls <- paste0(prefix, " = ", ifelse(sort(unique(vals)) > 0, "+", ""),
                 sprintf("%.2f", sort(unique(vals))))
  factor(labs, levels = lvls)
}


# ============================================================================
# SIMULATION 1
# ============================================================================
#
# dodge_var controls layout:
#   "tau"      (default): rows = accuracy, cols = prevalence, dodge = tau
#   "accuracy"          : rows = tau,      cols = prevalence, dodge = accuracy
#

plot_sim1 <- function(dir           = "simulation_data/sim1/",
                      y_col,
                      truth_col,
                      y_label       = "Value",
                      scale_y       = 1,
                      interval      = c("quantile", "sd"),
                      sd_se_filter  = NULL,
                      exclude_bad   = TRUE,
                      exclude_models = "M7",
                      model_order   = paste0("M", 1:6),
                      dodge_var     = c("tau", "accuracy"),
                      dodge_w       = 0.7,
                      base_size     = 12,
                      inc_legend    = TRUE,
                      inc_violin    = TRUE) {

  interval  <- match.arg(interval)
  dodge_var <- match.arg(dodge_var)
  y_col     <- rlang::enquo(y_col)
  truth_col <- rlang::enquo(truth_col)

  data <- load_sim(dir, exclude_bad, exclude_models, model_order)

  if (!is.null(sd_se_filter)) {
    data <- data %>% dplyr::filter(sd_se_logit == sd_se_filter)
  }

  d <- data %>%
    dplyr::mutate(
      acc_lab = paste0("Se/Sp = ", sprintf("%.2f", mean_se), "/", sprintf("%.2f", mean_sp)),
      prev_lab_raw = round(mu_prevalence * 100, 1),
      prev_lab     = paste0("\u03BC = ", sprintf("%.1f", prev_lab_raw)),
      prev_lab     = factor(prev_lab, levels = c("\u03BC = 1.5", "\u03BC = 6.0", "\u03BC = 12.0")),
      tau_lab      = paste0("\u03C4 = ", sprintf("%.1f", tau_logit)),
      tau_lab      = factor(tau_lab, levels = paste0("\u03C4 = ", sprintf("%.1f", sort(unique(tau_logit))))),
      tau_fac      = factor(sprintf("%.1f", tau_logit),
                            levels = sprintf("%.1f", sort(unique(tau_logit)))),
      y_val        = (!!y_col) * scale_y
    )

  if (dodge_var == "tau") {
    summ <- d %>%
      dplyr::group_by(acc_lab, prev_lab, model, tau_fac) %>%
      dplyr::summarise(
        mean_val = mean(y_val, na.rm = TRUE),
        sd_val   = stats::sd(y_val, na.rm = TRUE),
        q025     = stats::quantile(y_val, 0.025, na.rm = TRUE),
        q975     = stats::quantile(y_val, 0.975, na.rm = TRUE),
        truth    = dplyr::first((!!truth_col) * scale_y),
        .groups  = "drop"
      ) %>%
      dplyr::mutate(
        ymin  = if (interval == "sd") mean_val - sd_val else q025,
        ymax  = if (interval == "sd") mean_val + sd_val else q975,
        model = factor(model, levels = model_order)
      )

    hline_df <- summ %>% dplyr::distinct(acc_lab, prev_lab, truth)

    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                          linetype = "dashed", linewidth = 0.5, colour = "grey40")

    if (inc_violin) {
      p <- p + gghalves::geom_half_violin(
        data = d,
        ggplot2::aes(x = factor(model, levels = model_order), y = y_val,
                     group = interaction(model, tau_fac)),
        side = "r",
        trim = TRUE, scale = "width", alpha = 0.25,
        position = ggplot2::position_dodge(width = dodge_w)
      )
    }

    p <- p +
      ggplot2::geom_pointrange(
        data = summ,
        ggplot2::aes(x = model, y = mean_val, ymin = ymin, ymax = ymax,
                     shape = tau_fac, group = tau_fac),
        position = ggplot2::position_dodge(width = dodge_w), linewidth = 0.4
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(acc_lab), cols = ggplot2::vars(prev_lab)) +
      ggplot2::scale_shape_manual(name = expression(tau), values = c(17, 16))

  } else {
    # Build short accuracy labels for dodge
    d <- d %>%
      dplyr::mutate(
        acc_short = paste0(sprintf("%.2f", mean_se), "/", sprintf("%.2f", mean_sp)),
        acc_short = factor(acc_short, levels = unique(acc_short[order(mean_se)]))
      )

    summ <- d %>%
      dplyr::group_by(tau_lab, prev_lab, model, acc_short) %>%
      dplyr::summarise(
        mean_val = mean(y_val, na.rm = TRUE),
        sd_val   = stats::sd(y_val, na.rm = TRUE),
        q025     = stats::quantile(y_val, 0.025, na.rm = TRUE),
        q975     = stats::quantile(y_val, 0.975, na.rm = TRUE),
        truth    = dplyr::first((!!truth_col) * scale_y),
        .groups  = "drop"
      ) %>%
      dplyr::mutate(
        ymin  = if (interval == "sd") mean_val - sd_val else q025,
        ymax  = if (interval == "sd") mean_val + sd_val else q975,
        model = factor(model, levels = model_order)
      )

    hline_df <- summ %>% dplyr::distinct(tau_lab, prev_lab, truth)

    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                          linetype = "dashed", linewidth = 0.5, colour = "grey40")

    if (inc_violin) {
      p <- p + gghalves::geom_half_violin(
        data = d,
        ggplot2::aes(x = factor(model, levels = model_order), y = y_val,
                     group = interaction(model, acc_short)),
        side = "r",
        trim = TRUE, scale = "width", alpha = 0.25,
        position = ggplot2::position_dodge(width = dodge_w)
      )
    }

    p <- p +
      ggplot2::geom_pointrange(
        data = summ,
        ggplot2::aes(x = model, y = mean_val, ymin = ymin, ymax = ymax,
                     shape = acc_short, group = acc_short),
        position = ggplot2::position_dodge(width = dodge_w), linewidth = 0.4
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(tau_lab), cols = ggplot2::vars(prev_lab)) +
      ggplot2::scale_shape_manual(name = "Se/Sp", values = c(17, 16))
  }

  p <- p +
    ggplot2::labs(x = "Model", y = y_label) +
    sim_theme(base_size) +
    ggplot2::theme(
      legend.position      = if (inc_legend) c(0.5, 0.99) else "none",
      legend.justification = c("center", "top"),
      legend.direction     = "horizontal",
      legend.background    = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank()
    )

  p
}


# ============================================================================
# SIMULATION 2
# ============================================================================

plot_sim2 <- function(dir           = "simulation_data/sim2/",
                      y_col,
                      truth_col,
                      y_label       = "Value",
                      scale_y       = 1,
                      interval      = c("quantile", "sd"),
                      exclude_bad   = TRUE,
                      exclude_models = "M7",
                      model_order   = paste0("M", 1:6),
                      base_size     = 12,
                      inc_violin    = TRUE) {

  interval  <- match.arg(interval)
  y_col     <- rlang::enquo(y_col)
  truth_col <- rlang::enquo(truth_col)

  data <- load_sim(dir, exclude_bad, exclude_models, model_order)

  summ <- data %>%
    dplyr::group_by(sn, omega_se, model) %>%
    dplyr::summarise(
      mean_val = mean(!!y_col, na.rm = TRUE) * scale_y,
      sd_val   = stats::sd(!!y_col, na.rm = TRUE) * scale_y,
      q025     = stats::quantile(!!y_col, 0.025, na.rm = TRUE) * scale_y,
      q975     = stats::quantile(!!y_col, 0.975, na.rm = TRUE) * scale_y,
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      ymin  = if (interval == "sd") mean_val - sd_val else q025,
      ymax  = if (interval == "sd") mean_val + sd_val else q975,
      model = factor(model, levels = model_order)
    )

  hline_df <- data %>%
    dplyr::distinct(sn, omega_se, truth = !!truth_col) %>%
    dplyr::mutate(truth = truth * scale_y)

  lab_sn <- function(x) paste0("N = ", x)
  lab_omega <- function(x) paste0("\u03C9 = ", x)

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = model, y = mean_val)) +
    ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                        linetype = "dashed", linewidth = 0.5, colour = "grey40")

  if (inc_violin) {
    p <- p + gghalves::geom_half_violin(
      data = data %>% dplyr::mutate(mean_val = (!!y_col) * scale_y),
      ggplot2::aes(x = factor(model, levels = model_order), y = mean_val),
      side = "r",
      trim = TRUE, scale = "width", alpha = 0.25
    )
  }

  p <- p +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = ymin, ymax = ymax)) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(omega_se), cols = ggplot2::vars(sn),
      labeller = ggplot2::labeller(
        sn = ggplot2::as_labeller(lab_sn),
        omega_se = ggplot2::as_labeller(lab_omega)
      )
    ) +
    ggplot2::labs(x = "Model", y = y_label) +
    sim_theme(base_size)

  p
}


# ============================================================================
# SIMULATION 3
# ============================================================================

plot_sim3 <- function(dir           = "simulation_data/sim3/",
                      y_col,
                      truth_col,
                      y_label       = "Value",
                      scale_y       = 1,
                      interval      = c("quantile", "sd"),
                      exclude_bad   = TRUE,
                      exclude_models = "M7",
                      model_order   = paste0("M", 1:6),
                      base_size     = 12,
                      inc_violin    = TRUE) {

  interval  <- match.arg(interval)
  y_col     <- rlang::enquo(y_col)
  truth_col <- rlang::enquo(truth_col)

  data <- load_sim(dir, exclude_bad, exclude_models, model_order)

  summ <- data %>%
    dplyr::group_by(kappa_prior, sigma_bias, model) %>%
    dplyr::summarise(
      mean_val = mean(!!y_col, na.rm = TRUE) * scale_y,
      sd_val   = stats::sd(!!y_col, na.rm = TRUE) * scale_y,
      q025     = stats::quantile(!!y_col, 0.025, na.rm = TRUE) * scale_y,
      q975     = stats::quantile(!!y_col, 0.975, na.rm = TRUE) * scale_y,
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      ymin  = if (interval == "sd") mean_val - sd_val else q025,
      ymax  = if (interval == "sd") mean_val + sd_val else q975,
      model = factor(model, levels = model_order)
    )

  hline_df <- data %>%
    dplyr::distinct(kappa_prior, sigma_bias, truth = !!truth_col) %>%
    dplyr::mutate(truth = truth * scale_y)

  lab_kappa <- function(x) paste0("\u03BA = ", x)
  lab_sigma <- function(x) paste0("\u03C3_bias = ", x)

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = model, y = mean_val)) +
    ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                        linetype = "dashed", linewidth = 0.5, colour = "grey40")

  if (inc_violin) {
    p <- p + gghalves::geom_half_violin(
      data = data %>% dplyr::mutate(mean_val = (!!y_col) * scale_y),
      ggplot2::aes(x = factor(model, levels = model_order), y = mean_val),
      side = "r",
      trim = TRUE, scale = "width", alpha = 0.25
    )
  }

  p <- p +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = ymin, ymax = ymax)) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(kappa_prior), cols = ggplot2::vars(sigma_bias),
      labeller = ggplot2::labeller(
        kappa_prior = ggplot2::as_labeller(lab_kappa),
        sigma_bias  = ggplot2::as_labeller(lab_sigma)
      )
    ) +
    ggplot2::labs(x = "Model", y = y_label) +
    sim_theme(base_size)

  p
}


# ============================================================================
# SIMULATION 4a/4b: Systematic prior misspecification
# Faceted by delta_sp (cols, ordered -.1, 0, +.1) × mu_prevalence (rows)
# Dodged by delta_se (ordered -.1, 0, +.1)
# ============================================================================

plot_sim4 <- function(dir           = "simulation_data/sim4a/",
                      y_col,
                      truth_col,
                      y_label       = "Value",
                      scale_y       = 1,
                      interval      = c("quantile", "sd"),
                      exclude_bad   = TRUE,
                      exclude_models = NULL,
                      model_order   = c("M1", "M6"),
                      dodge_w       = 0.7,
                      base_size     = 12,
                      inc_violin    = TRUE,
                      inc_legend    = TRUE) {

  interval  <- match.arg(interval)
  y_col     <- rlang::enquo(y_col)
  truth_col <- rlang::enquo(truth_col)

  data <- load_sim(dir, exclude_bad, exclude_models, model_order)

  d <- data %>%
    dplyr::mutate(
      prev_lab = paste0("\u03BC = ", sprintf("%.1f", round(mu_prevalence * 100, 1))),
      prev_lab = factor(prev_lab, levels = c("\u03BC = 1.5", "\u03BC = 6.0", "\u03BC = 12.0")),
      dsp_lab  = ordered_delta_lab(delta_sp, prefix = "\u03B4Sp"),
      dse_fac  = {
        labs <- ifelse(delta_se > 0,
                       sprintf("+%.2f", delta_se),
                       sprintf("%.2f", delta_se))
        lvls <- ifelse(sort(unique(delta_se)) > 0,
                       sprintf("+%.2f", sort(unique(delta_se))),
                       sprintf("%.2f", sort(unique(delta_se))))
        factor(labs, levels = lvls)
      },
      y_val    = (!!y_col) * scale_y
    )

  summ <- d %>%
    dplyr::group_by(prev_lab, dsp_lab, model, dse_fac) %>%
    dplyr::summarise(
      mean_val = mean(y_val, na.rm = TRUE),
      sd_val   = stats::sd(y_val, na.rm = TRUE),
      q025     = stats::quantile(y_val, 0.025, na.rm = TRUE),
      q975     = stats::quantile(y_val, 0.975, na.rm = TRUE),
      truth    = dplyr::first((!!truth_col) * scale_y),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      ymin  = if (interval == "sd") mean_val - sd_val else q025,
      ymax  = if (interval == "sd") mean_val + sd_val else q975,
      model = factor(model, levels = model_order)
    )

  hline_df <- summ %>% dplyr::distinct(prev_lab, dsp_lab, truth)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                        linetype = "dashed", linewidth = 0.5, colour = "grey40")

  if (inc_violin) {
    p <- p + gghalves::geom_half_violin(
      data = d,
      ggplot2::aes(x = factor(model, levels = model_order), y = y_val,
                   group = interaction(model, dse_fac)),
      trim = TRUE, scale = "width", alpha = 0.25, side = "r",
      position = ggplot2::position_dodge(width = dodge_w)
    )
  }

  p <- p +
    ggplot2::geom_pointrange(
      data = summ,
      ggplot2::aes(x = model, y = mean_val, ymin = ymin, ymax = ymax,
                   shape = dse_fac, group = dse_fac),
      position = ggplot2::position_dodge(width = dodge_w), linewidth = 0.4
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(prev_lab), cols = ggplot2::vars(dsp_lab)) +
    ggplot2::scale_shape_manual(name = expression(delta[Se]),
                                values = c(17, 16, 15)) +
    ggplot2::labs(x = "Model", y = y_label) +
    sim_theme(base_size) +
    ggplot2::theme(
      legend.position      = if (inc_legend) c(0.5, 0.99) else "none",
      legend.justification = c("center", "top"),
      legend.direction     = "horizontal",
      legend.background    = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank()
    )

  p
}


# ============================================================================
# SIMULATION 5a/5b: Common prior centre
# Faceted by mean_se (rows) × mu_prevalence (cols)
# ============================================================================

plot_sim5 <- function(dir           = "simulation_data/sim5a/",
                      y_col,
                      truth_col,
                      y_label       = "Value",
                      scale_y       = 1,
                      interval      = c("quantile", "sd"),
                      exclude_bad   = TRUE,
                      exclude_models = NULL,
                      model_order   = c("M1", "M6"),
                      base_size     = 12,
                      inc_violin    = TRUE) {

  interval  <- match.arg(interval)
  y_col     <- rlang::enquo(y_col)
  truth_col <- rlang::enquo(truth_col)

  data <- load_sim(dir, exclude_bad, exclude_models, model_order)

  d <- data %>%
    dplyr::mutate(
      acc_lab  = paste0("Prior = ", sprintf("%.2f", mean_se), "/",
                        sprintf("%.2f", mean_se - 0.05)),
      prev_lab = paste0("\u03BC = ", sprintf("%.1f", round(mu_prevalence * 100, 1))),
      prev_lab = factor(prev_lab, levels = c("\u03BC = 1.5", "\u03BC = 6.0", "\u03BC = 12.0")),
      y_val    = (!!y_col) * scale_y
    )

  summ <- d %>%
    dplyr::group_by(acc_lab, prev_lab, model) %>%
    dplyr::summarise(
      mean_val = mean(y_val, na.rm = TRUE),
      sd_val   = stats::sd(y_val, na.rm = TRUE),
      q025     = stats::quantile(y_val, 0.025, na.rm = TRUE),
      q975     = stats::quantile(y_val, 0.975, na.rm = TRUE),
      truth    = dplyr::first((!!truth_col) * scale_y),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      ymin  = if (interval == "sd") mean_val - sd_val else q025,
      ymax  = if (interval == "sd") mean_val + sd_val else q975,
      model = factor(model, levels = model_order)
    )

  hline_df <- summ %>% dplyr::distinct(acc_lab, prev_lab, truth)

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = model, y = mean_val)) +
    ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                        linetype = "dashed", linewidth = 0.5, colour = "grey40")

  if (inc_violin) {
    p <- p + gghalves::geom_half_violin(
      data = d,
      ggplot2::aes(x = factor(model, levels = model_order), y = y_val),
      side = "r",
      trim = TRUE, scale = "width", alpha = 0.25
    )
  }

  p <- p +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = ymin, ymax = ymax)) +
    ggplot2::facet_grid(rows = ggplot2::vars(acc_lab), cols = ggplot2::vars(prev_lab)) +
    ggplot2::labs(x = "Model", y = y_label) +
    sim_theme(base_size)

  p
}


# ============================================================================
# SIMULATION 6: Gold-standard anchoring
# Faceted by mu_prevalence (rows) × delta_sp (cols, ordered -.1, 0, +.1)
# Dodged by prop_gold (symbols)
# ============================================================================

plot_sim6 <- function(dir           = "simulation_data/sim6/",
                      y_col,
                      truth_col,
                      y_label       = "Value",
                      scale_y       = 1,
                      interval      = c("quantile", "sd"),
                      exclude_bad   = TRUE,
                      exclude_models = NULL,
                      model_order   = c("M1", "M6"),
                      dodge_w       = 0.7,
                      base_size     = 12,
                      inc_violin    = TRUE,
                      inc_legend    = TRUE) {

  interval  <- match.arg(interval)
  y_col     <- rlang::enquo(y_col)
  truth_col <- rlang::enquo(truth_col)

  data <- load_sim(dir, exclude_bad, exclude_models, model_order)

  d <- data %>%
    dplyr::mutate(
      prev_lab = paste0("\u03BC = ", sprintf("%.1f", round(mu_prevalence * 100, 1))),
      prev_lab = factor(prev_lab, levels = c("\u03BC = 1.5", "\u03BC = 6.0", "\u03BC = 12.0")),
      dsp_lab  = ordered_delta_lab(delta_sp, prefix = "\u03B4Sp"),
      pgold_fac = factor(prop_gold,
                         labels = sprintf("%.0f%%", sort(unique(prop_gold)) * 100)),
      y_val    = (!!y_col) * scale_y
    )

  summ <- d %>%
    dplyr::group_by(prev_lab, dsp_lab, model, pgold_fac) %>%
    dplyr::summarise(
      mean_val = mean(y_val, na.rm = TRUE),
      sd_val   = stats::sd(y_val, na.rm = TRUE),
      q025     = stats::quantile(y_val, 0.025, na.rm = TRUE),
      q975     = stats::quantile(y_val, 0.975, na.rm = TRUE),
      truth    = dplyr::first((!!truth_col) * scale_y),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      ymin  = if (interval == "sd") mean_val - sd_val else q025,
      ymax  = if (interval == "sd") mean_val + sd_val else q975,
      model = factor(model, levels = model_order)
    )

  hline_df <- summ %>% dplyr::distinct(prev_lab, dsp_lab, truth)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                        linetype = "dashed", linewidth = 0.5, colour = "grey40")

  if (inc_violin) {
    p <- p + gghalves::geom_half_violin(
      data = d,
      ggplot2::aes(x = factor(model, levels = model_order), y = y_val,
                   group = interaction(model, pgold_fac)),
      side = "r",
      trim = TRUE, scale = "width", alpha = 0.25,
      position = ggplot2::position_dodge(width = dodge_w)
    )
  }

  p <- p +
    ggplot2::geom_pointrange(
      data = summ,
      ggplot2::aes(x = model, y = mean_val, ymin = ymin, ymax = ymax,
                   shape = pgold_fac, group = pgold_fac),
      position = ggplot2::position_dodge(width = dodge_w), linewidth = 0.4
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(prev_lab), cols = ggplot2::vars(dsp_lab)) +
    ggplot2::scale_shape_manual(name = expression(italic(p) * "(gold)"),
                                values = c(17, 16, 15)) +
    ggplot2::labs(x = "Model", y = y_label) +
    sim_theme(base_size) +
    ggplot2::theme(
      legend.position      = if (inc_legend) c(0.5, 0.99) else "none",
      legend.justification = c("center", "top"),
      legend.direction     = "horizontal",
      legend.background    = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank()
    )

  p
}


# ============================================================================
# SIMULATION 7: Moderator detection
# Faceted by slope_true (rows) × prior_sp (cols)
# Dodged by prop_gold (symbols)
# ============================================================================

plot_sim7_slope <- function(dir           = "simulation_data/sim7/",
                            y_col,
                            truth_col     = NULL,
                            y_label       = "Value",
                            interval      = c("quantile", "sd", "mc", "none"),
                            exclude_bad   = TRUE,
                            exclude_models = NULL,
                            model_order   = c("M1", "Gold", "AllGold", "M6"),
                            dodge_w       = 0.7,
                            base_size     = 12,
                            inc_violin    = TRUE,
                            inc_legend    = TRUE) {

  interval  <- match.arg(interval)
  y_col     <- rlang::enquo(y_col)

  truth_col <- rlang::enquo(truth_col)
  has_truth <- !rlang::quo_is_null(truth_col)
  if (!has_truth) truth_col <- rlang::quo(NA_real_)

  data <- load_sim(dir, exclude_bad, exclude_models, model_order)

  d <- data %>%
    dplyr::mutate(
      slope_lab  = paste0("\u03B2 = ", sprintf("%.1f", slope_true)),
      slope_lab  = factor(slope_lab,
                          levels = paste0("\u03B2 = ", sprintf("%.1f", sort(unique(slope_true))))),
      prior_lab  = paste0("Prior Sp = ", sprintf("%.2f", prior_sp)),
      pgold_fac  = factor(prop_gold,
                          labels = sprintf("%.0f%%", sort(unique(prop_gold)) * 100)),
      y_val      = !!y_col
    )

  summ <- d %>%
    dplyr::group_by(slope_lab, prior_lab, model, pgold_fac) %>%
    dplyr::summarise(
      mean_val = mean(y_val, na.rm = TRUE),
      sd_val   = stats::sd(y_val, na.rm = TRUE),
      q025     = stats::quantile(y_val, 0.025, na.rm = TRUE),
      q975     = stats::quantile(y_val, 0.975, na.rm = TRUE),
      n_val    = sum(!is.na(y_val)),
      truth    = dplyr::first(!!truth_col),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      # "mc": 95% Monte Carlo interval for a proportion such as power or Type I
      # error, mean +/- 1.96 * sqrt(p (1 - p) / n) with n the retained replicates
      mc_se = sqrt(mean_val * (1 - mean_val) / n_val),
      ymin  = switch(interval, sd = mean_val - sd_val, mc = mean_val - 1.96 * mc_se, none = mean_val, q025),
      ymax  = switch(interval, sd = mean_val + sd_val, mc = mean_val + 1.96 * mc_se, none = mean_val, q975),
      model = factor(model, levels = model_order)
    )

  hline_df <- summ %>% dplyr::distinct(slope_lab, prior_lab, truth)

  if (has_truth) {
    hline_df <- summ %>% dplyr::distinct(slope_lab, prior_lab, truth)
    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                          linetype = "dashed", linewidth = 0.5, colour = "grey40")
  } else {
    p <- ggplot2::ggplot()
  }

  if (inc_violin) {
    p <- p + gghalves::geom_half_violin(
      data = d,
      ggplot2::aes(x = factor(model, levels = model_order), y = y_val,
                   group = interaction(model, pgold_fac)),
      side = "r",
      trim = TRUE, scale = "width", alpha = 0.25,
      position = ggplot2::position_dodge(width = dodge_w)
    )
  }

  p <- p +
    ggplot2::geom_pointrange(
      data = summ,
      ggplot2::aes(x = model, y = mean_val, ymin = ymin, ymax = ymax,
                   shape = pgold_fac, group = pgold_fac),
      position = ggplot2::position_dodge(width = dodge_w), linewidth = 0.4
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(slope_lab), cols = ggplot2::vars(prior_lab)) +
    ggplot2::scale_shape_manual(name = expression(italic(p) * "(gold)"),
                                values = c(17, 16, 15)) +
    ggplot2::scale_x_discrete(labels = c(KnownSeSp = "Known")) +   # display label only
    ggplot2::labs(x = "Model", y = y_label) +
    sim_theme(base_size) +
    ggplot2::theme(
      legend.position      = if (inc_legend) c(0.25, 0.99) else "none",
      legend.justification = c("center", "top"),
      legend.direction     = "horizontal",
      legend.background    = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank()
    )

  p
}

# ---- Plotting function: Sim 8 prevalence ----------------------------------

plot_sim8_prevalence <- function(dir           = "simulation_data/sim8/",
                                 sn_filter     = NULL,
                                 interval      = c("quantile", "sd"),
                                 model_order   = c("KWGA", "Joint", "M6d", "M6j"),
                                 dodge_w       = 0.7,
                                 base_size     = 12,
                                 inc_violin    = FALSE,
                                 inc_legend    = TRUE,
                                 free_y        = TRUE) {

  interval <- match.arg(interval)

  results <- read_sim8_rds(dir = dir)

  if (is.null(results) || nrow(results) == 0) {
    stop("No Simulation 8 results found in ", dir)
  }

  results <- results %>% dplyr::filter(
    comp_maxrhat < 1.01, joint_maxrhat < 1.01,
    m6d_maxrhat < 1.01, m6j_maxrhat < 1.01,
    comp_ndiv < 1, joint_ndiv < 1,
    m6d_ndiv < 1, m6j_ndiv < 1
  )

  if (!is.null(sn_filter)) {
    results <- results %>% dplyr::filter(sn == sn_filter)
  }

  d <- results %>%
    dplyr::transmute(
      sn            = sn,
      mu_prevalence = mu_prevalence,
      mean_sp       = mean_sp,
      prop_gold     = prop_gold,
      scenario_id   = scenario_id,
      rep_id        = rep_id,
      KWGA  = kwga_prev_mean,
      Joint = joint_prev_mean,
      M6d   = m6d_prev_mean,
      M6j   = m6j_prev_mean
    ) %>%
    tidyr::pivot_longer(
      cols      = c(KWGA, Joint, M6d, M6j),
      names_to  = "model",
      values_to = "prev_est"
    ) %>%
    dplyr::mutate(
      prev_lab = paste0("\u03BC = ", sprintf("%.1f", round(mu_prevalence * 100, 1))),
      prev_lab = factor(prev_lab, levels = c("\u03BC = 1.5", "\u03BC = 6.0", "\u03BC = 12.0")),
      sp_lab   = paste0("Sp = ", sprintf("%.2f", mean_sp)),
      sp_lab   = factor(sp_lab, levels = paste0("Sp = ", sprintf("%.2f", sort(unique(mean_sp))))),
      pgold_fac = factor(prop_gold,
                         labels = sprintf("%.0f%%", sort(unique(prop_gold)) * 100)),
      model = factor(model, levels = model_order),
      y_val = prev_est * 100
    )

  summ <- d %>%
    dplyr::group_by(prev_lab, sp_lab, model, pgold_fac, mu_prevalence) %>%
    dplyr::summarise(
      mean_val = mean(y_val, na.rm = TRUE),
      sd_val   = stats::sd(y_val, na.rm = TRUE),
      q025     = stats::quantile(y_val, 0.025, na.rm = TRUE),
      q975     = stats::quantile(y_val, 0.975, na.rm = TRUE),
      truth    = dplyr::first(mu_prevalence) * 100,
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      ymin  = if (interval == "sd") mean_val - sd_val else q025,
      ymax  = if (interval == "sd") mean_val + sd_val else q975,
      model = factor(model, levels = model_order)
    )

  hline_df <- summ %>% dplyr::distinct(prev_lab, sp_lab, truth)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                        linetype = "dashed", linewidth = 0.5, colour = "grey40")

  if (inc_violin) {
    p <- p + gghalves::geom_half_violin(
      data = d,
      ggplot2::aes(x = model, y = y_val,
                   group = interaction(model, pgold_fac)),
      side = "r",
      trim = TRUE, scale = "width", alpha = 0.25,
      position = ggplot2::position_dodge(width = dodge_w)
    )
  }

  p <- p +
    ggplot2::geom_pointrange(
      data = summ,
      ggplot2::aes(x = model, y = mean_val, ymin = ymin, ymax = ymax,
                   shape = pgold_fac, group = pgold_fac),
      position = ggplot2::position_dodge(width = dodge_w), linewidth = 0.4
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(prev_lab),
                        cols = ggplot2::vars(sp_lab),
                        scales = if (free_y) "free_y" else "fixed") +
    ggplot2::scale_shape_manual(name = expression(italic(p) * "(gold)"),
                                values = c(17, 16, 15)) +
    ggplot2::labs(x = "Method", y = "Prevalence (%)") +
    sim_theme(base_size) +
    ggplot2::theme(
      legend.position      = if (inc_legend) c(0.5, 0.99) else "none",
      legend.justification = c("center", "top"),
      legend.direction     = "horizontal",
      legend.background    = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank()
    )

  if (!free_y) {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, 40))
  }

  p
}


# ---- Plotting function: Sim 8 Sp recovery --------------------------------

plot_sim8_sp_recovery <- function(dir          = "simulation_data/sim8/",
                                  sn_filter    = NULL,
                                  interval     = c("quantile", "sd"),
                                  dodge_w      = 0.7,
                                  base_size    = 12,
                                  inc_violin   = FALSE,
                                  inc_legend   = TRUE) {

  interval <- match.arg(interval)

  results <- read_sim8_rds(dir = dir)

  if (is.null(results) || nrow(results) == 0) {
    stop("No Simulation 8 results found in ", dir)
  }

  results <- results %>% dplyr::filter(
    comp_maxrhat < 1.01, joint_maxrhat < 1.01,
    m6d_maxrhat < 1.01, m6j_maxrhat < 1.01,
    comp_ndiv < 1, joint_ndiv < 1,
    m6d_ndiv < 1, m6j_ndiv < 1
  )

  if (!is.null(sn_filter)) {
    results <- results %>% dplyr::filter(sn == sn_filter)
  }

  d <- results %>%
    dplyr::transmute(
      sn          = sn,
      mean_sp     = mean_sp,
      mu_prevalence = mu_prevalence,
      prop_gold   = prop_gold,
      KWGA        = disc_weighted_sp,
      Joint       = joint_sp_mean
    ) %>%
    tidyr::pivot_longer(
      cols      = c(KWGA, Joint),
      names_to  = "method",
      values_to = "sp_est"
    ) %>%
    dplyr::mutate(
      prev_lab = paste0("\u03BC = ", sprintf("%.1f", round(mu_prevalence * 100, 1))),
      prev_lab = factor(prev_lab, levels = c("\u03BC = 1.5", "\u03BC = 6.0", "\u03BC = 12.0")),
      sp_lab   = paste0("True Sp = ", sprintf("%.2f", mean_sp)),
      sp_lab   = factor(sp_lab, levels = paste0("True Sp = ", sprintf("%.2f", sort(unique(mean_sp))))),
      pgold_fac = factor(prop_gold,
                         labels = sprintf("%.0f%%", sort(unique(prop_gold)) * 100)),
      method = factor(method, levels = c("KWGA", "Joint"))
    )

  summ <- d %>%
    dplyr::group_by(prev_lab, sp_lab, method, pgold_fac, mean_sp) %>%
    dplyr::summarise(
      mean_val = mean(sp_est, na.rm = TRUE),
      sd_val   = stats::sd(sp_est, na.rm = TRUE),
      q025     = stats::quantile(sp_est, 0.025, na.rm = TRUE),
      q975     = stats::quantile(sp_est, 0.975, na.rm = TRUE),
      truth    = dplyr::first(mean_sp),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      ymin = if (interval == "sd") mean_val - sd_val else q025,
      ymax = if (interval == "sd") mean_val + sd_val else q975
    )

  hline_df <- summ %>% dplyr::distinct(prev_lab, sp_lab, truth)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(data = hline_df, ggplot2::aes(yintercept = truth),
                        linetype = "dashed", linewidth = 0.5, colour = "grey40")

  if (inc_violin) {
    p <- p + gghalves::geom_half_violin(
      data = d,
      ggplot2::aes(x = method, y = sp_est,
                   group = interaction(method, pgold_fac)),
      side = "r",
      trim = TRUE, scale = "width", alpha = 0.25,
      position = ggplot2::position_dodge(width = dodge_w)
    )
  }

  p <- p +
    ggplot2::geom_pointrange(
      data = summ,
      ggplot2::aes(x = method, y = mean_val, ymin = ymin, ymax = ymax,
                   shape = pgold_fac, group = pgold_fac),
      position = ggplot2::position_dodge(width = dodge_w), linewidth = 0.4
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(prev_lab),
                        rows = ggplot2::vars(sp_lab)) +
    ggplot2::scale_shape_manual(name = expression(italic(p) * "(gold)"),
                                values = c(17, 16, 15)) +
    ggplot2::labs(x = "Method", y = "Estimated Sp") +
    sim_theme(base_size) +
    ggplot2::theme(
      legend.position      = if (inc_legend) c(0.5, 0.99) else "none",
      legend.justification = c("center", "top"),
      legend.direction     = "horizontal",
      legend.background    = ggplot2::element_blank(),
      legend.key           = ggplot2::element_blank()
    )

  p
}



# Make Plots --------------------------------------------------------------

## Simulation 1: Prevalence -----------------------
figure_s1a <- plot_sim1(
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100,
  sd_se_filter = 0.15, dodge_var = "tau"
) + ggtitle(expression("Simulation 1: " * sigma[Se/Sp] == 0.15)) + ylim(c(0, 40))

figure_s1b <- plot_sim1(
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100,
  sd_se_filter = 0.30, dodge_var = "tau"
) + ggtitle(expression("Simulation 1: " * sigma[Se/Sp] == 0.30)) + ylim(c(0, 40))

## Simulation 1: Tau -----------------------
figure_s4a <- plot_sim1(
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau), sd_se_filter = 0.15,
  dodge_var = "accuracy"
) + ggtitle(expression("Simulation 1 " * tau * ": " * sigma[Se/Sp] == 0.15))

figure_s4b <- plot_sim1(
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau), sd_se_filter = 0.30,
  dodge_var = "accuracy"
) + ggtitle(expression("Simulation 1 " * tau * ": " * sigma[Se/Sp] == 0.30))

## Simulation 2: Prevalence ---------------------------------
figure_s1c <- plot_sim2(
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100
) + ggtitle("Simulation 2") + ylim(c(0, 33))

## Simulation 2: Tau (dodge = accuracy, rows = tau) -----------------------
figure_s4c <- plot_sim2(
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau)
) + ggtitle("Simulation 2")

## Simulation 3: Prevalence -----------------------------
figure_s1d <- plot_sim3(
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100
) + ggtitle("Simulation 3") + ylim(c(0, 33))

## Simulation 3: Tau ----------------------------
figure_s4d <- plot_sim3(
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau)
) + ggtitle("Simulation 3")

## Simulation 4a: Prevalence ------------------------
figure_2a <- plot_sim4(
  dir = "simulation_data/sim4a/",
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100
) + ggtitle("Simulation 4a") + ylim(c(0, 40))

## Simulation 4a: Tau -------------------------
figure_s5a <- plot_sim4(
  dir = "simulation_data/sim4a/",
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau)
) + ggtitle("Simulation 4a")

## Simulation 4b: Prevalence ----------------------
figure_s2a <- plot_sim4(
  dir = "simulation_data/sim4b/",
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100
) + ggtitle("Simulation 4b") + ylim(c(0, 40))

## Simulation 4b: Tau ----------------------
figure_s6a <- plot_sim4(
  dir = "simulation_data/sim4b/",
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau)
) + ggtitle("Simulation 4b")

## Simulation 5a: Prevalence -------------------------
figure_s2b <- plot_sim5(
  dir = "simulation_data/sim5a/",
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100
) + ggtitle("Simulation 5a")

## Simulation 5a: Tau -------------------------
figure_s6b <- plot_sim5(
  dir = "simulation_data/sim5a/",
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau)
) + ggtitle("Simulation 5a")

## Simulation 5b: Prevalence ------------------------
figure_s2c <- plot_sim5(
  dir = "simulation_data/sim5b/",
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100
) + ggtitle("Simulation 5b")

## Simulation 5b: Tau ------------------------
figure_s6c <- plot_sim5(
  dir = "simulation_data/sim5b/",
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau)
) + ggtitle("Simulation 5b")

## Simulation 6: Prevalence ------------
figure_2b <- plot_sim6(
  y_col = prev_mean, truth_col = mu_prevalence,
  y_label = "Prevalence (%)", scale_y = 100
) + ggtitle("Simulation 6")

## Simulation 6: Tau ----------------
figure_s5b <- plot_sim6(
  y_col = tau_mean, truth_col = tau_logit,
  y_label = expression(tau)
) + ggtitle("Simulation 6")

## Simulation 7: Slope -----------------------
figure_3a <- plot_sim7_slope(
  dir          = "simulation_data/sim7/",
  y_col = slope_mean, truth_col = slope_true,
  model_order   = c("M1", "Gold", "KnownSeSp", "AllGold", "M6", "M8"),
  y_label = "Moderator slope (log-odds)"
)

figure_3b <- plot_sim7_slope(
  dir         = "simulation_data/sim7/",
  y_col       = slope_power,
  truth_col   = NULL,
  model_order = c("M1", "Gold", "KnownSeSp", "AllGold", "M6", "M8"),
  y_label     = "Power or Type I Error",
  inc_violin    = FALSE,
  interval      = "mc"   # 95% Monte Carlo intervals, sqrt(p (1 - p) / n)
) + ggplot2::coord_cartesian(ylim = c(0, 1))

## Simulation 8: Prevalence (N = 20) -------------------
figure_4a <- plot_sim8_prevalence(
  sn_filter = 20,
  interval  = "quantile",
  inc_violin = TRUE,
  free_y = FALSE
) + ggtitle("Simulation 8: Prevalence (N = 20)")

## Simulation 8: Prevalence (N = 40) -------------------
figure_4b <- plot_sim8_prevalence(
  sn_filter = 40,
  interval  = "quantile",
  inc_violin = TRUE,
  free_y = FALSE
) + ggtitle("Simulation 8: Prevalence (N = 40)")

## Simulation 8: Sp (N = 20) -------------------
figure_4c <- plot_sim8_sp_recovery(
  sn_filter = 20,
  inc_violin = TRUE
) + ggtitle("Simulation 8: Sp (N = 20)")

## Simulation 8: Sp (N = 40) -------------------
figure_4d <- plot_sim8_sp_recovery(
  sn_filter = 40,
  inc_violin = TRUE
) + ggtitle("Simulation 8: Sp (N = 40)")

# Main Text Figures ---------------------------

source("scripts/0. Setup/0.6 Submission Figure Export.R")  # journal-format TIFF copies (figures/submission/)

# Simulations 4a and 6 (Prevalence)
figure_2 = figure_2a / figure_2b
ggsave("figures/Figure_2.png", figure_2, width = 12, height = 12, dpi = 300)
save_submission_figure(figure_2, 2, width = 12, height = 12)

# Simulation 7 (Slope Recovery)
figure_3 <- figure_3a / figure_3b
ggsave("figures/Figure_3.png", figure_3, width = 8, height = 12, dpi = 300)
save_submission_figure(figure_3, 3, width = 8, height = 12)

# Simulation 8 (Prevalence)
figure_4 <- (figure_4a | figure_4b) / (figure_4c | figure_4d)
ggsave("figures/Figure_4.png", figure_4, width = 18, height = 12, dpi = 300)
save_submission_figure(figure_4, 4, width = 18, height = 12)

# Note: the case-study figures (manuscript Figures 5-7) are produced by 2.0 and
# saved under figures/mcma_case_study/

# Supplementary Figures -------------------------

# Simulations 1 - 3 (Prevalence)
figure_s1 <- (figure_s1a | figure_s1b) / (figure_s1c | figure_s1d)
ggsave("figures/Figure_S1.png", figure_s1, width = 18, height = 14, dpi = 300)
save_submission_figure(figure_s1, "S1", width = 18, height = 14)

# Simulations 4b and 5 (Prevalence)
figure_s2 = figure_s2a / (figure_s2b | figure_s2c)
ggsave("figures/Figure_S2.png", figure_s2, width = 10, height = 10, dpi = 300)
save_submission_figure(figure_s2, "S2", width = 10, height = 10)

# Note: Figure S3 is the Joint Posterior of the Three-level Case Study Model

# Simulations 1 - 3 (Tau)
figure_s4 <- (figure_s4a | figure_s4b) / (figure_s4c | figure_s4d)
ggsave("figures/Figure_S4.png", figure_s4, width = 18, height = 14, dpi = 300)
save_submission_figure(figure_s4, "S4", width = 18, height = 14)

# Simulations 4a and 6 (Tau)
figure_s5 <- figure_s5a / figure_s5b
ggsave("figures/Figure_S5.png", figure_s5, width = 12, height = 12, dpi = 300)
save_submission_figure(figure_s5, "S5", width = 12, height = 12)

# Simulations 4b and 5 (Tau)
figure_s6 = figure_s6a / (figure_s6b | figure_s6c)
ggsave("figures/Figure_S6.png", figure_s6, width = 10, height = 10, dpi = 300)
save_submission_figure(figure_s6, "S6", width = 10, height = 10)


# Main Text Table Summaries ---------------------------------------------------------

#
# Table 2
#

# Sim 4a: averaged across delta_se
tab_2a <- make_summary_table(
  dir          = "simulation_data/sim4a/",
  group_vars   = c("delta_sp", "mu_prevalence")
)

# Sim 6: averaged across prevalence
tab_2b <- make_summary_table(
  dir          = "simulation_data/sim6/",
  group_vars   = c("prop_gold", "delta_sp")
)

#
# Table 3
#

# Sim 7
tab_3 <- make_moderator_table(
  dir          = "simulation_data/sim7/",
  model_filter = c("M1", "Gold", "KnownSeSp", "AllGold", "M6", "M8"),
  group_vars   = c("slope_true", "prop_gold", "prior_sp")
) %>%
  dplyr::select(
    slope_true, prop_gold, prior_sp,
    M1_power, M1_coverage,
    Gold_power, Gold_coverage,
    KnownSeSp_power, KnownSeSp_coverage,
    AllGold_power, AllGold_coverage,
    M6_power, M6_coverage,
    M8_power, M8_coverage
  )

#
# Table 4
#

# Sim 8
tab_4 <- make_sim8_table(
  group_vars = c("mean_sp", "mu_prevalence")
) %>%
  dplyr::transmute(
    Sp_true       = mean_sp,
    mu_pi         = mu_prevalence,
    M8_prev_bias  = sprintf("%.3f", m8_prev_bias),
    M8_prev_cov   = sprintf("%.2f", m8_prev_coverage),
    M6j_prev_bias = sprintf("%.3f", m6j_prev_bias),
    M6j_prev_cov  = sprintf("%.2f", m6j_prev_coverage),
    M6d_prev_bias = sprintf("%.3f", m6d_prev_bias),
    M6d_prev_cov  = sprintf("%.2f", m6d_prev_coverage),
    KWGA_prev_bias = sprintf("%.3f", kwga_prev_bias),
    KWGA_prev_cov  = sprintf("%.2f", kwga_prev_coverage),
    M8_sp_bias    = sprintf("%.3f", joint_sp_bias),
    M8_sp_cov     = sprintf("%.2f", joint_sp_coverage),
    KWGA_sp_bias   = sprintf("%.3f", disc_sp_bias)
  ) %>%
  print(n = 20)

# Supplementary Table Summaries ---------------------------------------------------------

#
# Table S2
#

# Sim 1
tab_s2 <- make_summary_table(
  dir          = "simulation_data/sim1/",
  model_filter = c("M1", "M2", "M3", "M4", "M5", "M6"),
  group_vars   = c("mean_se", "mean_sp", "mu_prevalence", "tau_logit", "model")
) %>%
  tidyr::pivot_wider(
    names_from  = model,
    values_from = c(bias, coverage, ci_width, n),
    names_glue  = "{model}_{.value}"
  ) %>%
  dplyr::select(
    mean_se, mean_sp, mu_prevalence, tau_logit,
    M1_bias, M1_coverage,
    M2_bias, M2_coverage,
    M3_bias, M3_coverage,
    M4_bias, M4_coverage,
    M5_bias, M5_coverage,
    M6_bias, M6_coverage
  )

#
# Table S3
#

# Sim 2
tab_s3 <- make_summary_table(
  dir          = "simulation_data/sim2/",
  model_filter = c("M1", "M2", "M3", "M4", "M5", "M6"),
  group_vars   = c("sn", "omega_se", "model")
) %>%
  tidyr::pivot_wider(
    names_from  = model,
    values_from = c(bias, coverage, ci_width, n),
    names_glue  = "{model}_{.value}"
  ) %>%
  dplyr::select(
    sn, omega_se,
    M1_bias, M1_coverage,
    M2_bias, M2_coverage,
    M3_bias, M3_coverage,
    M4_bias, M4_coverage,
    M5_bias, M5_coverage,
    M6_bias, M6_coverage
  )

#
# Table S4
#

# Sim 3
tab_s4 <- make_summary_table(
  dir          = "simulation_data/sim3/",
  model_filter = c("M1", "M3", "M5", "M6"),
  group_vars   = c("sigma_bias", "kappa_prior", "model")
) %>%
  tidyr::pivot_wider(
    names_from  = model,
    values_from = c(bias, coverage, ci_width, n),
    names_glue  = "{model}_{.value}"
  ) %>%
  dplyr::select(
    sigma_bias, kappa_prior,
    M1_bias, M1_coverage,
    M3_bias, M3_coverage,
    M5_bias, M5_coverage,
    M6_bias, M6_coverage
  )

#
# Table S5
#

# Sim 4b
tab_s5a <- make_summary_table(
  dir          = "simulation_data/sim4b/",
  model_filter = "M6",
  group_vars   = c("delta_sp", "mu_prevalence")
)

# Sim 5a
tab_s5b <- make_summary_table(
  dir          = "simulation_data/sim5a/",
  model_filter = "M6",
  group_vars   = c("mean_se", "mean_sp", "mu_prevalence")
)

# Sim 5b
tab_s5c <- make_summary_table(
  dir          = "simulation_data/sim5b/",
  model_filter = "M6",
  group_vars   = c("mean_se", "mean_sp", "mu_prevalence")
)
