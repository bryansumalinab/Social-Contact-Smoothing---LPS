rm(list = ls())
source('cubicbs.R')
source('contactLPS.R')

data <- readRDS("data_holiday.RDS")


set.seed(0519)
model <- contactLPS(data = data,
                    K.age = 15,
                    K.disp = 10,
                    WAIC = TRUE)

model$WAIC

K.age <- model$K.age
K.disp <- model$K.disp

############# (A) Compute estimated contact rates ##############
############## for participant x contact age grid ##############
library(dplyr)
library(tidyr)
library(ggplot2)

### (A.1) Age grid
min_age_part <- min(data$age_part)
max_age_part <- max(data$age_part)
min_age_cont <- min(data$age_cont)
max_age_cont <- max(data$age_cont)
age_part_seq <- seq(min_age_part, max_age_part)
age_cont_seq <- seq(min_age_cont, max_age_cont)

grid_data <- expand.grid(
  age_part = age_part_seq,
  age_cont = age_cont_seq
)

### (A.2) B-spline basis matrices
build_symmetric_basis <- function(Bi, Bj, K) {
  Bgrid <- kronecker(Bj, Bi)
  Bmat <- matrix(1:(K^2), nrow = K, ncol = K)
  ind1 <- Bmat[lower.tri(Bmat)]
  ind2 <- t(Bmat)[lower.tri(t(Bmat))]
  Bgrid[,ind1] <- Bgrid[,ind1] + Bgrid[,ind2]
  Bsym <- Bgrid[,-c(Bmat[upper.tri(Bmat)])]
  Bsym
}

B_part <- cubicbs(age_part_seq, lower = min_age_part, upper = max_age_part, K = K.age)$Bmatrix
B_cont <- cubicbs(age_cont_seq, lower = min_age_cont, upper = max_age_cont, K = K.age)$Bmatrix
B_sym <- build_symmetric_basis(B_part, B_cont, K.age)

### (A.3) Model parameters
W_mat <- matrix(1, nrow = nrow(B_sym), ncol = 2)
W_groups <- model$W
group_names <- colnames(W_groups)
ximu_mode <- model$ximu_mode
pW <- model$pW
pB <- model$pB

xi_mu_W <- ximu_mode[seq_len(pW)]
xi_mu_B <- ximu_mode[-seq_len(pW)]
xi_B_list <- split(xi_mu_B, rep(seq_len(pW), each = pB))

### (A.4) Contact rate estimates
compute_rate <- function(group) {
  if (group == 1) {
    X <- cbind(W_mat[, 1], B_sym)
    coef <- c(xi_mu_W[1], xi_B_list[[1]])
  } else {
    X <- cbind(W_mat, B_sym, B_sym)
    coef <- c(xi_mu_W[c(1, group)], xi_B_list[[1]], xi_B_list[[group]])
  }
  exp(X %*% coef)
}

rates_list <- lapply(seq_along(xi_B_list), compute_rate)
rates_list <- setNames(rates_list, paste0("rate", seq_along(rates_list)))
rates_df <- cbind(grid_data, as.data.frame(rates_list))
rate_long <- rates_df %>%
  pivot_longer(
    cols = starts_with("rate"),
    names_to = "time",
    values_to = "rate"
  ) %>%
  mutate(
    time = as.integer(sub("rate", "", time))
  ) %>% 
  arrange(time, age_cont, age_part) %>% 
  mutate(group = factor(rep(group_names, each = length(age_part_seq)*length(age_cont_seq))))


### (A.5) Variance of log-contact rates
Covar_mu <- model$Covar[1:length(ximu_mode), 1:length(ximu_mode)]

compute_varlograte <- function(group) {
  if (group == 1) {
    idx <- c(1, seq(pW + 1, length.out = pB))
    X <- cbind(W_mat[, 1], B_sym)
    coef <- c(xi_mu_W[1], xi_B_list[[1]])
  } else {
    # indices for mu parameters
    mu_idx <- c(1, group, seq(pW + 1, length.out = pB), seq((group - 1)*pB + pW + 1, length.out = pB))
    X <- cbind(W_mat, B_sym, B_sym)
    coef <- c(xi_mu_W[c(1, group)], xi_B_list[[1]], xi_B_list[[group]])
    idx <- mu_idx
  }
  Cov_sub <- Covar_mu[idx, idx]
  var_log <- rowSums((X %*% Cov_sub) * X)
  mean_log <- as.numeric(X %*% coef)
  list(
    varlograte = var_log,
    varrate = (exp(var_log) - 1) * exp(2*mean_log + var_log)
  )
}


vars_list <- lapply(seq_along(xi_B_list), compute_varlograte)

varlog_list <- setNames(
  lapply(vars_list, `[[`, "varlograte"),
  paste0("varlograte", seq_along(vars_list))
)
varlog_df <- cbind(
  grid_data,
  as.data.frame(varlog_list)
)
varlog_long <- varlog_df %>%
  pivot_longer(
    cols = starts_with("varlograte"),
    names_to = "time",
    values_to = "varlograte"
  ) %>%
  mutate(time = as.integer(sub("varlograte", "", time))) %>%
  arrange(time, age_cont, age_part)

### (A.6) Estimated overdispersion
disp_est <- grid_data %>%
  mutate(
    disp     = as.numeric(model$disp_surface),
    inv_disp = 1/disp
  )

rate_long$group <- factor(rate_long$group,
                          levels = c("(Intercept)",
                                     "group1"))

################ (B) Plots ##################
############################################

gg_base <- theme_bw() +
  theme(
    axis.text      = element_text(size = 14),
    axis.title     = element_text(size = 14),
    legend.text    = element_text(size = 12),
    legend.title   = element_text(size = 14),
    strip.text     = element_text(size = 14, face = "bold")
  )

# Custom labels
group_labels <- c(`(Intercept)` = "No Holiday", group1 = "Holiday")

### (B.1) Contact rate surface
plot_lograte <- ggplot(rate_long, aes(age_part, age_cont, z = log(rate))) +
  geom_contour(color = "black") +
  geom_contour_filled(aes(fill = after_stat(level))) +
  scale_fill_viridis_d(option = "D", name = "log-rate") +
  facet_wrap(~ group, labeller = labeller(group = group_labels)) +
  labs(x = "Participant's age", y = "Contact's age") +
  gg_base

### (B.2) Overdispersion surface
plot_logdisp <- ggplot(disp_est, aes(age_part, age_cont, z = log(inv_disp))) +
  geom_contour() + geom_contour_filled() +
  scale_fill_viridis_d(name = "Overdispersion") +
  labs(x = "Participant's age", y = "Contact's age") + 
  theme_bw() +
  theme(
    axis.text      = element_text(size = 20),
    axis.title     = element_text(size = 20),
    legend.text    = element_text(size = 18),
    legend.title   = element_text(size = 20),
    strip.text     = element_text(size = 14, face = "bold")
  )

### (B.3) Mean difference and significance

combined_df <- rate_long %>%
  left_join(varlog_long, by = c("age_part", "age_cont", "time")) %>%
  mutate(meanlograte = as.numeric(log(rate)))
n_samples <- 1000
sample_lograte <- function(mean, var) {
  rnorm(n_samples, mean = mean, sd = sqrt(var))
}
lograte_samples <- t(apply(combined_df, 1, function(row) {
  sample_lograte(as.numeric(row["meanlograte"]), as.numeric(row["varlograte"]))
}))


block_size <- length(age_part_seq) * length(age_cont_seq)
n_groups   <- length(xi_B_list)

# row‑index vector for the intercept, repeated (n_groups-1) times
idx_int <- rep(1:block_size, times = n_groups - 1)

# row‑index vector for all non‑intercept rows
idx_grp <- seq_len(nrow(lograte_samples))[- (1:block_size)]

# Difference: group (holiday) - reference (no holiday)
diff_mat2 <- lograte_samples[idx_grp, , drop = FALSE] - 
  lograte_samples[idx_int, , drop = FALSE]


# mean diffreence, 2.5% and 97.5% quantiles
diff_all_summary <- t(apply(diff_mat2, 1, function(x) {
  c(mean = mean(x), quantile(x, 0.025), quantile(x, 0.975))
}))
diff_all_summary <- as.data.frame(diff_all_summary)

# Determine significance (whether 0 is outside the interval)
diff_all_summary$significant <- as.integer(diff_all_summary$`2.5%` > 0 | diff_all_summary$`97.5%` < 0)

diff_all <- cbind(rate_long[- (1:block_size),], diff_all_summary)


# Create a new column indicating direction and significance of difference
diff_all <- diff_all %>%
  mutate(rate_sign = ifelse(mean > 0, "Positive", "Negative")) %>% 
  mutate(rate_sign = factor(rate_sign,
                            levels = c("Positive",
                                       "Negative"),
                            labels = c("Positive",
                                       "Negative")
  )) %>% 
  mutate(sig_dir = case_when(
    !significant             ~ "Not significant",
    significant & mean > 0   ~ "Positive & significant",
    significant & mean < 0   ~ "Negative & significant"
  )) %>%
  mutate(sig_dir = factor(sig_dir,
                          levels = c("Positive & significant",
                                     "Negative & significant",
                                     "Not significant"),
                          labels = c("Yes (positive)",
                                     "Yes (negative)",
                                     "No")
  ))


plot_meandiff <- ggplot(diff_all, aes(age_part, age_cont, z = mean)) +
  geom_contour(color = "black") + geom_contour_filled(aes(fill = after_stat(level))) +
  scale_fill_viridis_d(option = "D", name = "Mean\ndifference") +
  labs(x = "Participant's age", y = "Contact's age") + 
  theme_bw() +
  theme(
    axis.text      = element_text(size = 20),
    axis.title     = element_text(size = 20),
    legend.text    = element_text(size = 20),
    legend.title   = element_text(size = 20),
    strip.text     = element_text(size = 14, face = "bold")
  )


plot_direction <- ggplot(diff_all, aes(age_part, age_cont, z = mean)) +
  geom_contour(color = "black") + geom_contour_filled(aes(fill = rate_sign)) +
  scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato"),
                    name = "Direction of\nthe difference") +
  labs(x = "Participant's age", y = "Contact's age") + 
  theme_bw() +
  theme(
    axis.text      = element_text(size = 20),
    axis.title     = element_text(size = 20),
    legend.text    = element_text(size = 20),
    legend.title   = element_text(size = 20),
    strip.text     = element_text(size = 14, face = "bold")
  )



plot_significance <- ggplot(diff_all, aes(age_part, age_cont, fill = sig_dir)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    name   = "Significance of\nthe difference\nat 5% level",
    values = c(
      "No"             = "grey80",
      "Yes (positive)" = "#4472c4",
      "Yes (negative)"   = "tomato"
    )
  ) +
  labs(x = "Participant's age", y = "Contact's age") + 
  theme_bw() +
  theme(
    axis.text      = element_text(size = 20),
    axis.title     = element_text(size = 20),
    legend.text    = element_text(size = 20),
    legend.title   = element_text(size = 20),
    strip.text     = element_text(size = 14, face = "bold")
  )



plot_lograte # contact rate surface
plot_logdisp # overdispersion surface
plot_meandiff # Differences in contact rates (Holiday - no Holiday)
plot_direction # Direction of the difference
plot_significance # Significance of the difference

