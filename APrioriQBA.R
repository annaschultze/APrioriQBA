
# SCRIPT INFORMATION 
# ------------------------------------------------
# script: a_priori_QBA.R
# author: A Schultze
# description: QBA on power calculations
#
# ------------------------------------------------
# 0. HOUSEKEEPING --------------------------------

library(tidyverse)
library(gt)
library(ggplot2)

# 1. FUNCTIONS --------------------------------

#' A priori QBA outcome misclassification
#' 
#' Conduct simple QBA for a priori outcome misclassification 
#' Note: this uses an input HR (taken as true) and recalculates on an OR (for ease of formula application)
#' Assumption is that these will be close as long as outcome is reasonable rare 
#' 
#' @param N Total sample size - does not matter
#' @param p_exp Prevalence of exposure - does not matter
#' @param rate-unexp Incidence rate in the unexposed 
#' @param fu Amount of FU (for final outcome numbers)
#' @param hr Observed (target) HR - assumed to be true
#' @param se_exp1 sensitivity in the exposed 
#' @param se_exp0 sensitivity in the unexposed 
#' @param sp_exp1 specificity in the exposed 
#' @param sp_exp0 specificity in the unexposed 
#' 
#' @return df with true and OR that would be observed (obs_OR), by rate

apriori_qba_outcome_misclassification <- function(N, 
                                                 p_exp, 
                                                 rate_unexp,
                                                 fu,
                                                 hr, 
                                                 se_exp1, 
                                                 se_exp0, 
                                                 sp_exp1, 
                                                 sp_exp0) { 
  
  # as per Fox, events in rows, exposure in columns
  a_raw <- N * p_exp * rate_unexp * fu * hr # exposed cases
  c_raw <- N * p_exp - a_raw # exposed non-cases 
  b_raw <- N * (1 - p_exp) * rate_unexp *fu # unexposed cases
  d_raw <- N * (1 - p_exp) - b_raw # unexposed non-cases 
  
  # 'true' RR 
  raw_OR <- (a_raw / c_raw) / (b_raw / d_raw)
  
  a_corr <- a_raw * se_exp1 + c_raw * (1 - sp_exp1)
  b_corr <- b_raw * se_exp0 + d_raw * (1 - sp_exp0)
  c_corr <- a_raw * (1 - se_exp1) + c_raw * sp_exp1
  d_corr <- b_raw * (1 - se_exp0) + d_raw * sp_exp0
  
  if (any(c(a_corr, b_corr, c_corr, d_corr) < 0 )) { 
    
    warning("negative cell counts")
    
    obs_OR <- "negative cell counts"
    
  } else {
    
    obs_OR <- (a_corr / b_corr) / (c_corr / d_corr)
    obs_OR <- round(obs_OR, 3)
    obs_OR <- as.character(obs_OR)
    
  }

  data.frame(fu   = fu,
             se_exp1      = se_exp1, 
             se_exp0      = se_exp0,
             sp_exp1      = sp_exp1, 
             sp_exp0      = sp_exp0,
             true_OR       = round(raw_OR,  3),
             obs_OR      = obs_OR
  )

}

#' A priori QBA exposure misclassification
#' 
#' 'Reverse' QBA to calculate observed, from an assumed truth  
#' Conduct simple QBA for exposure misclassification 
#' Note: this uses an input HR (assumed true) and recalculates on an OR (for ease of formula application)
#' Assumption is that these will be close as long as outcome is reasonable rare 
#' 
#' @param N Total sample size - does not matter, but can vary to confirm this
#' @param p Prevalence of exposure - does not matter, but can vary to confirm this
#' @param rate-unexp Incidence rate in the unexposed (impacts total number of cases)
#' @param hr True HR (assumed)
#' @param fu Amount of follow-up (impacts total number of cases)
#' @param se_outc1 sensitivity in cases 
#' @param se_outc0 sensitivity in non-cases 
#' @param sp_outc1 specificity in cases 
#' @param sp_outc0 specificity in non-cases
#' 
#' @return df with true and OR that would be observed (obs_OR), by rate

apriori_qba_exposure_misclassification <- function(N, 
                                                 p_exp,
                                                 fu, 
                                                 rate_unexp, 
                                                 hr, 
                                                 se_outc1, 
                                                 se_outc0, 
                                                 sp_outc1, 
                                                 sp_outc0) { 
  
  # as per Fox, events in rows, exposure in columns
  a_raw <- N * p_exp * rate_unexp * fu * hr # exposed cases
  c_raw <- N * p_exp - a_raw # exposed non-cases 
  b_raw <- N * (1 - p_exp) * rate_unexp *fu # unexposed cases
  d_raw <- N * (1 - p_exp) - b_raw # unexposed non-cases 
  
  # 'true' RR 
  raw_OR <- (a_raw / c_raw) / (b_raw / d_raw)
  
  a_corr <- a_raw * se_outc1 + b_raw * (1 - sp_outc1)
  b_corr <- a_raw * (1 - se_outc1) + b_raw * sp_outc1
  c_corr <- c_raw * se_outc0 + d_raw * (1 - sp_outc0)
  d_corr <- c_raw * (1 - se_outc0) + d_raw * sp_outc0
  
  if (any(c(a_corr, b_corr, c_corr, d_corr) < 0 )) { 
    
    warning("negative cell counts")
    
    obs_OR <- "negative cell counts"
    
  } else {
    
    obs_OR <- (a_corr / c_corr) / (b_corr / d_corr)
    obs_OR <- round(obs_OR, 3)
    obs_OR <- as.character(obs_OR)
    
  }
  
  data.frame(fu   = fu,
             se_outc1      = se_outc1, 
             se_outc0      = se_outc0,
             sp_outc1      = sp_outc1, 
             sp_outc0      = sp_outc0,
             true_OR       = round(raw_OR,  3),
             obs_OR      = obs_OR
  )
  
}

#' A priori QBA unmeasured confounding
#' 
#' Conduct a priori QBA for unmeasured confounding
#' Note: this uses an input HR but formula is for RR 
#' Assumption is that these will be close as long as outcome is reasonable rare 
#' 
#' @param p_exp1 Prevalence of confounder among exposed
#' @param rr_ux Ratio of confounder prevalence in the unexposed relative to the exposed (p_exp0 / p_exp1)
#' @param rr_uyx rate ratio for association between confounder and outcome
#' @param hr assumed true (target) HR
#' 
#' @return df with true and the RR that would be observed, obs_RR, by prevalence


apriori_qba_unmeasured_confounding <- function(p_exp1, 
                                               rr_ux,
                                               rr_uy, 
                                               hr) { 
  
  p_exp0 <- rr_ux * p_exp1 
  
  qba_numerator <- (p_exp1 * (rr_uy - 1)) + 1
  qba_denominator <- (p_exp0 * (rr_uy - 1)) + 1
  
  obs_rr <- hr * (qba_numerator / qba_denominator)
  
  data.frame(p_exp1   = p_exp1,
             p_exp0   = p_exp0, 
             rr_ux    = rr_ux, 
             rr_uy    = rr_uy,
             hr       = hr,
             obs_rr  = round(obs_rr,  3))
  
}

# 2. APPLY TO SCENARIOS --------------------------------

# outcome misclassification scenarios 
N <- 10000
p_exp <- 0.30
hr = 1.2
rate_unexp = 0.015

fu_vals <- c(3, 5, 10)
sens_vals <- c(0.4, 0.5, 0.6, 0.7, 0.75, 0.80)
spec_vals <- c(0.75, 0.80, 0.85, 0.90, 0.95, 0.99)

sens_df <- expand_grid(
  fu = fu_vals,
  sens = sens_vals) %>%
  mutate(
    N = N, 
    p_exp = p_exp, 
    hr = hr,
    rate_unexp = rate_unexp, 
    se_exp1 = sens, 
    se_exp0 = sens,
    sp_exp1 = 1,
    sp_exp0 = 1) %>% 
  select(-sens)

spec_df <- expand_grid(
  fu = fu_vals,
  spec = spec_vals) %>%
  mutate(
    N = N, 
    p_exp = p_exp, 
    hr = hr,
    rate_unexp = rate_unexp, 
    se_exp1 = 1, 
    se_exp0 = 1,
    sp_exp1 = spec,
    sp_exp0 = spec) %>% 
  select(-spec)

results_outc_sens <- purrr::pmap_dfr(sens_df, apriori_qba_outcome_misclassification)
results_outc_spec <- purrr::pmap_dfr(spec_df, apriori_qba_outcome_misclassification)

write.csv(results_outc_sens, "results_outc_sens.csv", row.names = FALSE)
write.csv(results_outc_spec, "results_outc_spec.csv", row.names = FALSE)

# unmeasured confounding 
hr <- 1.2

p_exp1_vals  <- c(0.15, 0.20, 0.25)
rr_ux_vals <- c(1.15, 1.20, 1.25, 1.50, 2.00)
rr_uy_vals <- c(1.2, 1.5, 2.0, 4.0, 6.0)

uc_frailty_df <- expand_grid(
  p_exp1 = p_exp1_vals,
  rr_ux = rr_ux_vals, 
  rr_uy = rr_uy_vals) %>% 
  mutate(hr = hr)

results_uc_frailty <- purrr::pmap_dfr(uc_frailty_df, apriori_qba_unmeasured_confounding)
write.csv(results_uc_frailty, "results_uc_frailty.csv", row.names = FALSE)

# exposure misclassification scenarios 
N <- 10000
p_exp <- 0.3
hr = 1.2
rate_unexp = 0.015

fu_vals <- c(3, 5, 10)
spec_vals <- c(0.5, 0.6, 0.7, 0.8, 0.9)

spec_df <- expand_grid(
  fu = fu_vals,
  spec = spec_vals) %>%
  mutate(
    N = N, 
    p_exp = p_exp, 
    hr = hr,
    rate_unexp = rate_unexp, 
    se_outc1 = 1, 
    se_outc0 = 1,
    sp_outc1 = spec,
    sp_outc0 = spec) %>% 
  select(-spec)

results_exp_spec <- purrr::pmap_dfr(spec_df, apriori_qba_exposure_misclassification)

write.csv(results_exp_spec, "results_exp_spec.csv", row.names = FALSE)

# exposure misclassification scenarios - differential
N <- 10000
p_exp <- 0.3
hr = 1.2
rate_unexp = 0.015

fu_vals <- c(3, 5, 10)
spec_vals <- c(0.5, 0.6, 0.7, 0.8, 0.9)

spec_df <- expand_grid(
  fu = fu_vals,
  spec = spec_vals) %>%
  mutate(
    N = N, 
    p_exp = p_exp, 
    hr = hr,
    rate_unexp = rate_unexp, 
    se_outc1 = 1, 
    se_outc0 = 1,
    sp_outc1 = spec,
    sp_outc0 = spec+0.1) %>% 
  select(-spec)

results_exp_spec_diff <- purrr::pmap_dfr(spec_df, apriori_qba_exposure_misclassification)

write.csv(results_exp_spec_diff, "results_exp_spec_diff.csv", row.names = FALSE)

# 4. PLOTS --------------------------------------------

results_outc_sens$baseline_label <- factor(paste0("Years of follow-up = ", results_outc_sens$fu),
                            levels = paste0("Years of follow-up = ", c(3, 5, 10)))

results_outc_sens$obs_OR <- as.numeric(results_outc_sens$obs_OR)

ggplot(results_outc_sens) +
  geom_point(aes(x = se_exp1, y = obs_OR), colour = "skyblue3") +
  geom_line(aes(x = se_exp1, y = true_OR), colour = "lavenderblush3", linetype = "dashed") +
  facet_wrap(~ baseline_label) +
  labs(
    x = "Sensitivity",
    y = "Observed Odds Ratio",
    title = "Outcome Mislcassification: Impact of low sensitivity",
    caption = "Dashed line = True OR"
  ) +
  scale_y_continuous(limits = c(1.150, 1.3500)) +
  theme(panel.spacing.x = unit(2, "lines"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "gray95"))

ggsave("outc_sensitivity.png", width = 10, height = 5, dpi = 300)


# Specificity 
results_outc_spec$invalid <- results_outc_spec$obs_OR == "negative cell counts"
results_outc_spec$obs_OR_num <- as.numeric(ifelse(results_outc_spec$invalid, NA, results_outc_spec$obs_OR))
results_outc_spec$true_OR <- as.numeric(results_outc_spec$true_OR)

results_outc_spec$baseline_label <- factor(paste0("Years of follow-up = ", results_outc_spec$fu),
                                      levels = paste0("Years of follow-up = ", c(3, 5, 10)))

results_outc_spec$sp_exp1 <- factor(results_outc_spec$sp_exp1)

# Data frame of invalid specificity values for shading
invalid_sp <- results_outc_spec[results_outc_spec$invalid, c("baseline_label", "sp_exp1", "true_OR")]

ggplot(results_outc_spec, aes(x = sp_exp1)) +
  geom_point(aes(y = obs_OR_num), colour = "skyblue3") +
  geom_hline(data = unique(results_outc_spec[, c("baseline_label", "true_OR")]),
             aes(yintercept = true_OR), colour = "lavenderblush3", linetype = "dashed") +
  geom_point(data = invalid_sp, aes(x = sp_exp1, y = true_OR),
             shape = 4, size = 3, colour = "grey40") +
  facet_wrap(~ baseline_label, scales = "free_y") +
  labs(
    x = "Specificity",
    y = "Observed Odds Ratio",
    title = "Outcome Misclassification: Impact of Low Specificity",
    caption = "Dashed line = True OR; × = negative cell counts (invalid)"
  ) +
  scale_y_continuous(limits = c(0.800, 1.500), labels = scales::label_number(accuracy = 0.01)) +
  theme(panel.spacing.x = unit(2, "lines"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "gray95"))

ggsave("outc_specificity.png", width = 10, height = 5, dpi = 300)

# Unmeasured Confounding - Frailty 
results_uc_frailty$rr_ux <- factor(results_uc_frailty$rr_ux)
results_uc_frailty$rr_uy <- factor(results_uc_frailty$rr_uy)

ggplot(results_uc_frailty, aes(x = rr_uy, y = obs_rr, colour = rr_ux, group = rr_ux)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1.2, linetype = "dashed", colour = "lavenderblush3") +
  facet_wrap(~ p_exp1, labeller = labeller(p_exp1 = function(x) paste0("Prevalence in exposed = ", x))) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01), limits = c(0.7, 1.3)) +
  labs(
    x = "RR: Frailty - Dementia",
    y = "Corrected RR",
    colour = "RR: Prevalence of frailty in unexposed compared to exposed",
    title = "QBA: Unmeasured Confounding",
    caption = "Dashed line = True HR"
  ) +
  theme(panel.spacing.x = unit(2, "lines"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "gray95"), 
        legend.position = "bottom") + 
  scale_colour_brewer(palette = "BrBG")

ggsave("uc_frailty.png", width = 10, height = 5, dpi = 300)

# Exposure misclassification - non differential 

results_exp_spec$baseline_label <- factor(paste0("Years of follow-up = ", results_exp_spec$fu),
                                           levels = paste0("Years of follow-up = ", c(3, 5, 10)))

results_exp_spec$obs_OR <- as.numeric(results_exp_spec$obs_OR)

ggplot(results_exp_spec) +
  geom_point(aes(x = sp_outc1, y = obs_OR), colour = "skyblue3") +
  geom_line(aes(x = sp_outc1, y = true_OR), colour = "lavenderblush3", linetype = "dashed") +
  facet_wrap(~ baseline_label) +
  labs(
    x = "Specificity",
    y = "Observed Odds Ratio",
    title = "Exposure Mislcassification: Impact of low specificity",
    caption = "Dashed line = True OR"
  ) +
  scale_y_continuous(limits = c(1.00, 2.00)) +
  theme(panel.spacing.x = unit(2, "lines"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "gray95"))

ggsave("exp_specificity.png", width = 10, height = 5, dpi = 300)

# Exposure misclassification - differential 

results_exp_spec_diff$baseline_label <- factor(paste0("Years of follow-up = ", results_exp_spec_diff$fu),
                                          levels = paste0("Years of follow-up = ", c(3, 5, 10)))

results_exp_spec_diff$obs_OR <- as.numeric(results_exp_spec_diff$obs_OR)

ggplot(results_exp_spec_diff) +
  geom_point(aes(x = sp_outc1, y = obs_OR), colour = "skyblue3") +
  geom_line(aes(x = sp_outc1, y = true_OR), colour = "lavenderblush3", linetype = "dashed") +
  facet_wrap(~ baseline_label) +
  labs(
    x = "Specificity",
    y = "Observed Odds Ratio",
    title = "Exposure Mislcassification: Impact of low specificity",
    caption = "Dashed line = True OR"
  ) +
  scale_y_continuous(limits = c(1.00, 2.00)) +
  theme(panel.spacing.x = unit(2, "lines"),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "gray95"))

ggsave("exp_specificity_diff.png", width = 10, height = 5, dpi = 300)




