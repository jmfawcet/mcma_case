# NOTE:
#
# The observed prevalence under misclassification is:
#   p_obs = Se * θ + (1 - Sp)(1 - θ)
#
# Bias is the difference between what is observed and what is true:
#   bias = p_obs - θ
# = Se * θ + (1 - Sp)(1 - θ) - θ
# = Se * θ + 1 - Sp - θ + Sp * θ - θ
# = (Se + Sp - 2)θ + (1 - Sp)
# = (1 - Sp) - (1 - Sp)θ - (1 - Se)θ
# = (1 - Sp)(1 - θ) - (1 - Se)θ
# where (1 - Sp) is FPR and (1 - Se) is FNR,
# and (1 - θ), θ are the truly-negative and truly-positive fractions:
# = FPR * (1 - θ) - FNR * θ

# Bias is zero when (i.e., the cross over point):
#
#   θ = FPR / (FPR + FNR)
#
# Below this, false positives dominate and prevalence is overestimated. Above this,
# false negatives dominate and it's underestimated. For typical psychiatric screeners
# (Se = .85, Sp = .80), θ0 = .20/(.20 + .15) = .57 well above any realistic disorder
# prevalence, so overestimation is the norm.


# Load libraries ----------------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)


# Linear Association (Not in Paper)  --------------------------------

prev_seq = seq(0.001, 0.999, by = 0.001)

configs = tibble(
  label = c("Se=.95, Sp=.95", "Se=.90, Sp=.85",
            "Se=.85, Sp=.80", "Se=.85, Sp=.75",
            "Se=.80, Sp=.70", "Se=.75, Sp=.65"),
  se = c(.95, .90, .85, .85, .80, .75),
  sp = c(.95, .85, .80, .75, .70, .65)
)

bias_df = expand_grid(true_prev = prev_seq, config = seq_len(nrow(configs))) %>%
  mutate(
    se    = configs$se[config],
    sp    = configs$sp[config],
    label = configs$label[config],
    fpr   = 1 - sp,
    fnr   = 1 - se,
    bias  = fpr * (1 - true_prev) - fnr * true_prev,
    crossover = fpr / (fpr + fnr)
  )

crossover_df = bias_df %>%
  distinct(label, crossover, se, sp)

ggplot(bias_df, aes(x = true_prev, y = bias * 100, colour = label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(data = crossover_df,
             aes(x = crossover, y = 0),
             size = 2.5, shape = 16) +
  scale_x_continuous(
    name = "True Prevalence",
    labels = scales::percent_format(),
    breaks = seq(0, 1, by = 0.10)
  ) +
  labs(
    y = "Bias (difference)",
    colour = "Diagnostic Accuracy"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

# Heat Map (Figure 1)  --------------------------------

prev_seq = seq(0.01, 0.99, by = 0.01)
sp_seq   = seq(0.50, 1.00, by = 0.01)
se_vals  = c(.5, .6, 0.70, 0.80, 0.90, 1.00)

bias_df = expand_grid(
  true_prev = prev_seq,
  sp        = sp_seq,
  se        = se_vals
) %>%
  mutate(
    fpr      = 1 - sp,
    fnr      = 1 - se,
    bias     = fpr * (1 - true_prev) - fnr * true_prev,
    se_label = sprintf("Se = %.2f", se),
    se_label = factor(se_label, levels = sprintf("Se = %.2f", se_vals))
  )

fig1 = ggplot(bias_df, aes(x = true_prev, y = sp, fill = bias * 100)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(data = bias_df %>% filter(se < 1),
               aes(x = true_prev, y = sp, z = bias * 100),
               inherit.aes = FALSE,
               breaks = 0,
               colour = "black", linewidth = 0.8, linetype = "dashed") +
  # scale_fill_gradient2(
  #   low = "#2166AC", mid = "white", high = "#B2182B",
  #   midpoint = 0, name = "Bias"
  # ) +
  # scale_fill_gradientn(
  #   colours = c("#2166AC", "#67A9CF", "#D1E5F0", "white",
  #               "#FDDBC7", "#EF8A62", "#B2182B"),
  #   values  = scales::rescale(c(-1, -0.25, -0.05, 0, 0.05, 0.25, 1)),
  #   name    = "Bias",
  #   limits  = c(-50, 50)
  # ) +
  scale_fill_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "white",
                "#FDDBC7", "#EF8A62", "#B2182B"),
    values  = scales::rescale(c(-1, -0.25, -0.05, 0, 0.05, 0.25, 1)),
    name    = "Bias",
    limits  = c(-20, 20),
    oob     = scales::squish
  ) +
  scale_x_continuous(
    name = "True Prevalence",
    labels = scales::percent_format(),
    breaks = seq(0, 1, by = 0.20)
  ) +
  scale_y_continuous(
    name = "Specificity (Sp)",
    breaks = seq(0.50, 1.00, by = 0.10)
  ) +
  facet_wrap(~ se_label, ncol = 3) +

  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

#ggsave('figures/Figure_1.pdf', plot=fig1, width=12, height=6)
ggsave('figures/Figure_1.png', plot=fig1, width=12, height=6)
source("scripts/0. Setup/0.6 Submission Figure Export.R")  # journal-format TIFF copies (figures/submission/)
save_submission_figure(fig1, 1, width = 12, height = 6)
