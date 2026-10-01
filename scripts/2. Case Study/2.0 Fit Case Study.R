####
# Case Study: Major Depressive Disorder in Public Safety Personnel
# =================================================================
#
# Companion analysis for the `mcma` package (Bayesian
# misclassification-corrected meta-analysis).
#
# Pipeline:
#   1. Naive (uncorrected) prevalence benchmark [M1]
#   2. Gold-standard-only anchor
#   3. Comparison model + discrete KWGA over Se/Sp
#   4. Joint model (M8) estimating prevalence, Se, and Sp
#      - Two-level (study) and three-level (study within report)
#   5. Corrected model (M6) with priors centred at M8 estimates
#   6. Sensitivity grid of M6 fits across Se x Sp prior centres
#   7. Posterior predictive checks
#
# Data: White, N., Wagner, S. L., Matthews, L. R., Randall, C., Regehr, C., White, M., Alden, L. E., Buys, N., 
# Carey, M. G., Corneil, W., Fyfe, T., Fraess-Phillips, A., & Krutop, E. (2025). Methodological Correlates of 
# Variability in Depressive and Anxiety Symptoms in Public Safety Personnel: A Systematic Review and 
# Metaregression. Traumatology, 31(2), 346–361. https://doi.org/10.1037/trm0000538
# 
# Parameterization: bounded Se/Sp (Se, Sp in [0.5, 1)).
#
# NOTE: The sensitivity grid takes a really long time; results are
# cached under `models/case_study/`. Delete that directory to refit from
# scratch and set RERUN = TRUE.
####


# Setup --------------------------------------------------------------------

# remotes::install_github("jmfawcet/mcma")
library(mcma) # Note this may overwrite some functions used in the earlier simulation code
library(brms)
library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(patchwork)

options(brms.backend = "cmdstanr")
set.seed(12318238) # For reproducibility

fig_dir   <- "figures/mcma_case_study"
model_dir <- "models/mcma_case_study"
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function that runs a specific expression if it hasn't 
# already been run and saved to disk; rerun makes it ALWAYS
# rerun the expression.
cache <- function(name, expr, rerun = FALSE) {
  fp   <- file.path(model_dir, paste0(name, ".rds"))
  call <- substitute(expr)
  if (is.call(call) && is.name(call[[1]]) &&
      startsWith(as.character(call[[1]]), "mcma_fit")) {
    call$file       <- fp
    call$file_refit <- if (rerun) "always" else "on_change"
    return(eval(call, parent.frame()))
  }
  if (!rerun && file.exists(fp)) return(readRDS(fp))
  out <- expr
  saveRDS(out, fp)
  out
}

PREV_CENTER <- 0.06  # prior centre on population MDD prevalence
RERUN <- FALSE # Set to TRUE to rerun models even if cached
ITER = 20000 # Number of iterations
WARMUP = 5000 # Number of warmups

# Data ---------------------------------------------------------------------

# Download the case-study data from OSF (White et al., 2025)
# OSF project: https://osf.io/9y2xj/  (depression in public safety personnel)
data_file <- "data/dep_psp.csv"
if (!file.exists(data_file)) {
  dir.create("data", recursive = TRUE, showWarnings = FALSE)
  download.file("https://osf.io/download/9y2xj/", data_file,
                mode = "wb", quiet = TRUE)
}

# Read in and process the data file
dat <- read.csv("data/dep_psp.csv") %>%
  mutate(
    y     = round(DepPrev * analyticsamplesize),
    n     = analyticsamplesize,
    es_id = factor(row_number()),
    DepressionTool = recode( # In one case the original data coded PHQ9 as PHQ, which we believe was a typo, fixing
      DepressionTool,
      PHQ = "PHQ9"
    ), 
    cutoff_type = case_when(
      Depthreshold == "DSMcrit"                            ~ "diagnostic",
      Depthreshold == "gpnorm"                             ~ "norm_referenced",
      Depthreshold == "notreported" | is.na(Depthreshold)  ~ "unknown",
      TRUE                                                 ~ "self_report_numeric"
    ),
    is_gold        = DepressionTool2 == "diagnostic",   
    measure_cutoff = factor(paste(DepressionTool, Depthreshold, sep = "_"))
  ) %>%
  filter(catquality != "low")

measure_map <- dat %>%
  filter(!is_gold) %>%
  distinct(measure_cutoff) %>%
  arrange(measure_cutoff) %>%
  mutate(measure_id = as.character(row_number()))

dat <- dat %>%
  left_join(measure_map, by = "measure_cutoff") %>%
  mutate(measure_id = ifelse(is_gold, "gold", measure_id))

cat(sprintf(
  "Estimates: %d  |  Gold: %d  |  Screen: %d  |  Distinct measures: %d | Distinct measure + cutoff combinations: %d\n",
  nrow(dat), sum(dat$is_gold), sum(!dat$is_gold), n_distinct(dat$DepressionTool), n_distinct(dat$measure_cutoff)
))


# Priors -------------------------------------------------------------------

# Uniform Se = Sp = 0.80 with effective sample size 100 for all screening
# measures; gold-standard measures are fixed to near-perfect accuracy. These
# centres serve as the base priors for M6 and are overwritten per grid point
# in the sensitivity analysis.
screen_measures <- sort(unique(dat$measure_id[!dat$is_gold]))
all_measures    <- c(screen_measures, "gold")

base_priors <- mcma_priors(
  measure_id = all_measures,
  se         = c(rep(0.80, length(screen_measures)), NA),
  sp         = c(rep(0.80, length(screen_measures)), NA),
  kappa      = c(rep(100,  length(screen_measures)), NA),
  is_gold    = c(rep(FALSE, length(screen_measures)), TRUE)
)


# M1: Naive benchmark ------------------------------------------------------

fit_naive <- cache("M1_naive_case",
  mcma_fit(
    data        = dat,
    priors      = NULL,
    prev_center = PREV_CENTER,
    prev_re     = ~ (1 | es_id),
    iter = ITER,
    warmup = WARMUP
  ), RERUN
)


# Gold-standard-only anchor ------------------------------------------------

fit_gold <- cache("gold_only_case",
  mcma_fit(
    data        = filter(dat, is_gold),
    priors      = NULL,
    prev_center = PREV_CENTER,
    prev_re     = ~ (1 | es_id),
    iter = ITER,
    warmup = WARMUP
  ), RERUN
)


# Comparison model + discrete KWGA over Se x Sp -----------------------------

# Centre of the prior for the screeners
scr_rate <- with(dat, sum(y[!is_gold]) / sum(n[!is_gold]))

fit_comp <- cache("comparison_gold_screen_case",
  mcma_fit_comparison(
    data        = dat,
    gold_column = "is_gold",
    gold_prev_center = PREV_CENTER,
    scr_prev_center  = scr_rate,
    prev_re     = ~ (1 | es_id),
    iter = ITER,
    warmup = WARMUP
  ), RERUN
)

kwga      <- mcma_kwga(fit_comp,
                     se_grid = seq(0.70, .975, by = 0.025),
                     sp_grid = seq(0.70, .975, by = 0.025),
                     clamp_scoring = FALSE)
kwga_prev <- mcma_kwga_prevalence(kwga, fit_comparison = fit_comp, n_draws = 4000,
                                  resample = "joint", seed = 12318238)

print(kwga)
cat(sprintf("KWGA prevalence: %.2f%% [%.2f%%, %.2f%%]\n",
            kwga_prev$summary$mean  * 100,
            kwga_prev$summary$ci_lb * 100,
            kwga_prev$summary$ci_ub * 100))

# M8: Joint model (two-level) ----------------------------------------------

fit_m8 <- cache("M8_joint_case",
  mcma_fit_joint(
    data            = dat,
    gold_column     = "is_gold",
    prev_center     = PREV_CENTER,
    se_prior_center = 0.80, sp_prior_center = 0.80,
    se_prior_sd     = 1.0,  sp_prior_sd     = 1.0,
    prev_re         = ~ (1 | es_id),
    sesp_re         = ~ (1 | es_id),
    correlated_re   = TRUE,
    bounded         = TRUE,
    iter = ITER,
    warmup = WARMUP
  ), RERUN
)

m8_prev <- extract_prevalence(fit_m8)
m8_sesp <- extract_sesp(fit_m8, level = "global", bounded=TRUE)
cat(sprintf(
  "M8: prev = %.2f%% [%.2f%%, %.2f%%]  Se = %.3f [%.3f, %.3f]  Sp = %.3f [%.3f, %.3f]\n",
  m8_prev$mean  * 100, m8_prev$l95  * 100, m8_prev$u95  * 100,
  m8_sesp$se_mean, m8_sesp$se_ci_lb, m8_sesp$se_ci_ub,
  m8_sesp$sp_mean, m8_sesp$sp_ci_lb, m8_sesp$sp_ci_ub
))

source("scripts/0. Setup/0.6 Submission Figure Export.R")  # journal-format TIFF copies (figures/submission/)

fig_7 <- mcma::plot_joint_sesp(fit_m8, n_draws = 20000, show_marginals = TRUE, bounded=TRUE, plot_title=NULL)
ggsave(file.path(fig_dir, "Figure_7.png"), fig_7, width = 6, height = 4)
save_submission_figure(fig_7, 7, width = 6, height = 4)


# M8: Joint model (three-level: studies within reports) --------------------

fit_m8_tl <- cache("M8_joint_threelevel_case",
  mcma_fit_joint(
    data          = dat,
    gold_column   = "is_gold",
    prev_center   = PREV_CENTER,
    prev_re       = ~ (1 | es_id) + (1 | refid),
    sesp_re       = ~ (1 | es_id) + (1 | refid),
    correlated_re = TRUE,
    bounded       = TRUE,
    iter = ITER,
    warmup = WARMUP
  ), RERUN
)

m8_tl_prev <- extract_prevalence(fit_m8_tl)
m8_tl_sesp <- extract_sesp(fit_m8_tl, level = "global", bounded=TRUE)
cat(sprintf(
  "M8 (3-level): prev = %.2f%% [%.2f%%, %.2f%%]  Se = %.3f [%.3f, %.3f]  Sp = %.3f [%.3f, %.3f]\n",
  m8_tl_prev$mean  * 100, m8_tl_prev$l95  * 100, m8_tl_prev$u95  * 100,
  m8_tl_sesp$se_mean, m8_tl_sesp$se_ci_lb, m8_tl_sesp$se_ci_ub,
  m8_tl_sesp$sp_mean, m8_tl_sesp$sp_ci_lb, m8_tl_sesp$sp_ci_ub
))

fig_s3 <- mcma::plot_joint_sesp(fit_m8_tl, n_draws = 20000, show_marginals = TRUE, bounded=TRUE, plot_title=NULL)
ggsave(file.path(fig_dir, "Figure_S3.png"), fig_s3, width = 6, height = 4)
save_submission_figure(fig_s3, "S3", width = 6, height = 4)

# M6: Corrected model with M8-informed priors ------------------------------

m6_priors_from_m8 <- update(base_priors,
                            se = m8_sesp$se_mean,
                            sp = m8_sesp$sp_mean)

fit_m6_m8 <- cache("M6_from_M8_case",
                   mcma_fit(
                     data          = dat,
                     priors        = m6_priors_from_m8,
                     prev_center   = PREV_CENTER,
                     prev_re       = ~ (1 | es_id),
                     sesp_re       = ~ (1 | es_id),
                     correlated_re = TRUE,
                     bounded       = TRUE,
                     iter = ITER,
                     warmup = WARMUP
                   ), RERUN
)

# Sensitivity grid ---------------------------------------------------------

se_vals <- seq(0.70, .95, by = 0.05)
sp_vals <- seq(0.70, .95, by = 0.05)

sens <- mcma_sensitivity(
  data          = dat,
  priors        = base_priors,
  se_grid       = se_vals,
  sp_grid       = sp_vals,
  prev_center   = PREV_CENTER,
  bounded       = TRUE,
  prev_re       = ~ (1 | es_id),
  sesp_re       = ~ (1 | es_id),
  correlated_re = TRUE,
  chains        = 4, cores = 4,
  iter = ITER, 
  warmup = WARMUP,
  model_dir     = file.path(model_dir, "sensitivity_case"),
  rerun = RERUN
)

sens_comp <- mcma_sensitivity_comparison(
  data          = dat,
  priors        = base_priors,
  se_grid       = se_vals,
  sp_grid       = sp_vals,
  prev_center   = PREV_CENTER,
  bounded       = TRUE,
  prev_re       = ~ (1 | es_id),
  sesp_re       = ~ (1 | es_id),
  correlated_re = TRUE,
  gold_column   = "is_gold",
  model_dir     = file.path(model_dir, "sensitivity_comp_case"),
  chains        = 4, cores = 4,
  iter          = ITER,
  warmup = WARMUP,
  rerun = RERUN
)

# For highlighting
sp_levels <- levels(factor(sens$results$prior_sp))
se_levels <- levels(factor(sens$results$prior_se))

g_prev <- plot(sens,      type = "heatmap",    metric = "prev_mean", text_size = 2.2) +
  ggplot2::annotate("rect",
                    xmin = match("0.8", sp_levels) - 0.5,
                    xmax = match("0.9", sp_levels) + 0.5,
                    ymin = match("0.8", se_levels) - 0.5,
                    ymax = match("0.9", se_levels) + 0.5,
                    fill = NA, colour = "red", linewidth = 1.5
                    )

g_diff <- plot(sens_comp, type = "difference", include_cell_estimates = TRUE, text_size = 2.0) +
  ggplot2::annotate("rect",
                    xmin = match("0.8", sp_levels) - 0.5,
                    xmax = match("0.9", sp_levels) + 0.5,
                    ymin = match("0.8", se_levels) - 0.5,
                    ymax = match("0.9", se_levels) + 0.5,
                    fill = NA, colour = "red", linewidth = 1.5
  )

fig_6 <- (g_prev / g_diff)
ggsave(file.path(fig_dir, "Figure_6.png"), fig_6, width = 12, height = 12)
save_submission_figure(fig_6, 6, width = 12, height = 12)

# Forest plot of prevalence estimates --------------------------------------

prev_table <- bind_rows(
  extract_prevalence(fit_naive) %>% mutate(model = "M1 (Naive)"),
  extract_prevalence(fit_gold)  %>% mutate(model = "Gold-only"),
  extract_prevalence(fit_m6_m8) %>% mutate(model = "M6 (M8-informed)"),
  extract_prevalence(fit_m8)    %>% mutate(model = "M8 (Joint)"),
  tibble(
    model  = "KWGA (analytic)",
    mean   = kwga_prev$summary$mean,
    median = kwga_prev$summary$median,
    sd     = kwga_prev$summary$sd,
    l95    = kwga_prev$summary$ci_lb,
    u95    = kwga_prev$summary$ci_ub
  )
) %>% dplyr::select(model, mean, median, l95, u95)

grid_subset <- sens$results %>%
  filter(prior_se %in% seq(0.70, .95, by = 0.05),
         prior_sp %in% seq(0.70, .95, by = 0.05)) %>%
  transmute(
    model  = sprintf("M6 (Se=%.2f, Sp=%.2f)", prior_se, prior_sp),
    mean   = prev_mean,
    median = prev_median,
    l95    = ci_lb,
    u95    = ci_ub
  )

highlighted <- c("M1 (Naive)", "Gold-only", "M8 (Joint)",
                 "M6 (M8-informed)", "KWGA (analytic)")

forest_df <- bind_rows(prev_table, grid_subset) %>%
  mutate(model = factor(model, levels = model[order(mean)]))

forest_plot <- ggplot(forest_df, aes(x = mean * 100, y = model)) +
  geom_pointrange(aes(xmin = l95 * 100, xmax = u95 * 100),
                  linewidth = 0.4, fatten = 2) +
  geom_point(data = ~ filter(.x, model %in% highlighted),
             size = 4, shape = 18, colour = "firebrick") +
  scale_x_continuous("Prevalence (%)", breaks = seq(0, 20, 2)) +
  labs(y = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

ggsave(file.path(fig_dir, "Figure_5.png"),
       forest_plot, width = 7, height = 6)
save_submission_figure(forest_plot, 5, width = 7, height = 6)

cat("\n--- Prevalence Estimates ---\n")
print(prev_table %>%
        mutate(across(c(mean, median, l95, u95),
                      ~ sprintf("%.2f%%", . * 100))), n = 20)


# Convergence diagnostics --------------------------------------------------

conv <- mcma_convergence(
  M1         = fit_naive,
  Gold_only  = fit_gold,
  Comparison = fit_comp,
  M8         = fit_m8,
  M8_3level  = fit_m8_tl,
  M6_M8      = fit_m6_m8
)
print(conv)


# Posterior predictive checks: primary models -----------------------------

ppc_table <- tibble(
  model = c("M1 (Naive)", "Gold-only", "M8 (Joint)", "M6 (M8-informed)"),
  fit   = list(fit_naive, fit_gold, fit_m8, fit_m6_m8)
) %>%
  mutate(ppc = map(fit, mcma_ppc)) %>%
  unnest(ppc) %>%
  dplyr::select(-fit, -n_obs)

print(ppc_table)

