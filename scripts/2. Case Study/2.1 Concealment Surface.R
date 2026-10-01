# 2.1 Concealment Surface ---------------------------------------------------
#
# Partial identification surface for the concealment extension of M8 (see the
# Discussion). c is the probability that a true case (successfully) conceals the 
# disorder in a diagnostic interview, so interview sensitivity becomes (1 - c); 
# gamma is the fraction of that concealment carried over to screening measures, 
# whose sensitivity becomes (1 - c * gamma) times its baseline. Neither is
# identified by a single estimate per study, so both are fixed and M8 is refit
# over a grid of assumed values; c = 0 is M8 itself. The corrected prevalence
# from each fit traces theta_hat(c, gamma), drawn as Figure S7.
# 
# Needs to run 2.0 first, uses stuff from from that script. Takes a *very* long time to 
# run!

if (!exists("dat") || !exists("cache")) {
  stop("Run 2.0 Fit Case Study.R first; 2.1 uses the data and helpers it creates.")
}

SURFACE_ITER     <- ITER   # iterations per chain, half warm-up, as in the main case-study fits
SURFACE_C        <- seq(0, 0.9, by = 0.1)
SURFACE_GAMMA    <- c(0, 0.25, 0.50, 0.75, 1)
SURFACE_SEED     <- 12318238
RUN_SURFACE_FITS <- TRUE   # FALSE summarizes the cached cells only

surface_dir <- file.path(model_dir, "concealment_surface")
dir.create(surface_dir, recursive = TRUE, showWarnings = FALSE)


# Cell fitter ---------------------------------------------------------------

# Fit one (c, gamma) cell with the M8 specification from 2.0 and keep its
# prevalence summary and convergence diagnostics. The seed is fixed per cell
# and the chain length is recorded with the summary.
fit_conceal_cell <- function(c, gamma, iter = SURFACE_ITER, rerun = RERUN) {
  name <- sprintf("concealment_surface/c%03d_g%03d", round(100 * c), round(100 * gamma))
  cache(name, {
    fit <- mcma_fit_joint(
      data            = dat,
      gold_column     = "is_gold",
      prev_center     = PREV_CENTER,
      se_prior_center = 0.80, sp_prior_center = 0.80,
      se_prior_sd     = 1.0,  sp_prior_sd     = 1.0,
      prev_re         = ~ (1 | es_id),
      sesp_re         = ~ (1 | es_id),
      correlated_re   = TRUE,
      bounded         = TRUE,
      c               = c,
      gamma           = gamma,
      iter            = iter,
      warmup          = iter / 2,
      seed            = SURFACE_SEED + round(100 * c) + 1000 * round(100 * gamma),
      refresh         = 0
    )
    prev <- extract_prevalence(fit)
    conv <- mcma_convergence(fit = fit)
    tibble(
      c = c, gamma = gamma,
      mean = prev$mean, median = prev$median, l95 = prev$l95, u95 = prev$u95,
      max_rhat = conv$max_rhat, n_divergences = conv$n_divergences,
      min_ess_bulk = conv$min_ess_bulk, min_ess_tail = conv$min_ess_tail,
      iter = iter
    )
  }, rerun)
}

# gamma is inert when c = 0, so that cell is fitted once and reused below.
surface_cells <- expand.grid(c = SURFACE_C, gamma = SURFACE_GAMMA) %>%
  filter(c > 0 | gamma == SURFACE_GAMMA[1]) %>%
  arrange(gamma, c)

if (RUN_SURFACE_FITS) {
  for (i in seq_len(nrow(surface_cells))) {
    fit_conceal_cell(surface_cells$c[i], surface_cells$gamma[i])
  }
}


# Convergence repair --------------------------------------------------------

# Cells that fail the convergence gate used throughout the paper (R-hat > 1.01
# or any divergence) are refit once with chains twice as long; the longer fit
# replaces the cached summary. Cells that still fail are reported below rather
# than dropped.
SURFACE_ITER_REPAIR <- 2 * SURFACE_ITER

if (RUN_SURFACE_FITS) {
  for (iter_repair in SURFACE_ITER_REPAIR) {
    poor <- list.files(surface_dir, pattern = "\\.rds$", full.names = TRUE) %>%
      map_dfr(readRDS) %>%
      filter(max_rhat > 1.01 | n_divergences > 0, iter < iter_repair)
    for (i in seq_len(nrow(poor))) {
      fit_conceal_cell(poor$c[i], poor$gamma[i], iter = iter_repair, rerun = TRUE)
    }
  }
}


# Surface summary -----------------------------------------------------------

surface_fits <- list.files(surface_dir, pattern = "\\.rds$", full.names = TRUE) %>%
  map_dfr(readRDS)

# add c = 0
surface <- bind_rows(
  filter(surface_fits, c > 0),
  crossing(select(filter(surface_fits, c == 0), -gamma), gamma = SURFACE_GAMMA)
) %>%
  arrange(gamma, c)

# Reference estimates from earlier scripts
naive_prev     <- prev_table$mean[prev_table$model == "M1 (Naive)"]
gold_only_prev <- prev_table$mean[prev_table$model == "Gold-only"]

cat("\n--- Concealment surface: corrected prevalence (%) by assumed c (rows) and gamma (columns) ---\n")
surface %>%
  transmute(c, gamma, value = sprintf("%.2f", 100 * mean)) %>%
  pivot_wider(names_from = gamma, values_from = value, names_prefix = "gamma = ") %>%
  print(n = Inf)

cat(sprintf("\nReference estimates: M8 %.1f%%, gold-only model %.1f%%, naive model %.1f%%\n",
            100 * m8_prev$mean, 100 * gold_only_prev, 100 * naive_prev))
cat(sprintf("Convergence: %d of %d fitted cells with R-hat > 1.01 or divergences; max R-hat %.3f; min bulk ESS %.0f\n",
            sum(surface_fits$max_rhat > 1.01 | surface_fits$n_divergences > 0), nrow(surface_fits),
            max(surface_fits$max_rhat), min(surface_fits$min_ess_bulk)))


# Tipping points ------------------------------------------------------------

# The assumed c where the corrected estimate first reaches a target, by
# linear interpolation between the two neighbouring grid cells (NA when the
# grid never reaches it).
tipping_c <- function(surf, target) {
  surf <- arrange(surf, c)
  i <- which(surf$mean >= target)[1]
  if (is.na(i)) return(NA_real_)
  if (i == 1) return(surf$c[1])
  approx(surf$mean[(i - 1):i], surf$c[(i - 1):i], xout = target)$y
}

tipping <- surface %>%
  group_by(gamma) %>%
  group_modify(~ tibble(c_naive = tipping_c(.x, naive_prev))) %>%
  ungroup()

cat("\n--- Tipping point: assumed c at which the corrected estimate reaches the naive estimate ---\n")
print(tipping, n = Inf)
cat("(NA = not reached within the fitted grid.)\n")


# Figure S7 -----------------------------------------------------------------

# Corrected prevalence against assumed concealment, one panel per gamma: the
# grey ribbon is the 95% credible interval and the black line the posterior
# mean. The dashed line is M8 without concealment and the dotted line the naive
# estimate, labelled inside every panel.
source("scripts/0. Setup/0.6 Submission Figure Export.R")

# The naive label sits above its line; the corrected label sits just below the
# lowest edge of the ribbon, so it never overlaps the interval.
surface_labels <- crossing(
  gamma_lab = factor(sprintf("γ = %.2f", SURFACE_GAMMA)),
  tibble(y     = c(naive_prev, min(surface$l95) - 0.002),
         label = c("Naive Estimate", "Corrected Estimate"),
         vjust = c(-0.96, 1))
)

fig_s7 <- surface %>%
  mutate(gamma_lab = factor(sprintf("γ = %.2f", gamma))) %>%
  ggplot(aes(c)) +
  geom_ribbon(aes(ymin = l95, ymax = u95), fill = "grey72") +
  geom_hline(yintercept = m8_prev$mean, linetype = "dashed", colour = "grey15", linewidth = 0.4) +
  geom_hline(yintercept = naive_prev, linetype = "dotted", colour = "grey15", linewidth = 0.55) +
  geom_line(aes(y = mean), linewidth = 0.45) +
  geom_point(aes(y = mean), size = 0.7) +
  geom_text(data = surface_labels, aes(x = -0.02, y = y, label = label, vjust = vjust),
            hjust = 0, size = 5.5 / .pt, fontface = "bold") +
  facet_wrap(~ gamma_lab, nrow = 1) +
  scale_x_continuous("Assumed Concealment (c)",
                     breaks = c(0, .3, .6, .9), labels = c("0", ".3", ".6", ".9")) +
  scale_y_continuous("Corrected Prevalence",
                     breaks = seq(0, .30, .05), labels = paste0(seq(0, 30, 5), "%")) +
  coord_cartesian(ylim = c(-0.022, max(0.30, max(surface$u95)))) +
  theme_classic(base_size = 8) +
  theme(panel.border    = element_rect(fill = NA, colour = "black", linewidth = 0.5),
        axis.line       = element_blank(),
        strip.background = element_blank(),
        strip.text      = element_text(face = "bold", size = 8),
        panel.spacing.x = unit(3, "pt"),
        axis.text       = element_text(size = 7, colour = "black"))

ggsave(file.path(fig_dir, "Figure_S7.png"), fig_s7, width = 6.5, height = 1.9)
save_submission_figure(fig_s7, "S7", width = 6.5, height = 1.9)
