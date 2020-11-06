library(rjags) # to synthesise and simulate
library(R2jags) # to synthesise and simulate
library(dplyr) # to manipulate data
library(reshape2) # to transform data
library(ggplot2) # for plots
library(scales) # for dollar signs and commas
library(dampack) # for CEA and to calculate ICERs
library(tidyverse) # for general data wrangling and vis
library(reshape2) # data wrangling
library(BCEA) # bayesian CEA package
library(bayesplot) # tools for posterior inspection

# ==========================================================================================
# Evidence Synthesis ------------------------------------------------------
# ==========================================================================================
# Data on the probability of reponse (sub-model 1) ------------------------
# Number of studies
N_res <- c(2, 4)
# Cases for t = 0
r_res1 <- c(55, 42)
# Total n observations for t = 1
n_res1 <- c(165, 118)

# Cases for t = 1
r_res2 <- c(77, 61, 61, 25)
# Total n observations for t = 2
n_res2 <- c(161, 203, 143, 46)

# Data on probability of toxicity (sub-model 2) ---------------------------
# Number of studies
N_tox <- c(2, 4)
# Cases for t = 0
r_tox1 <- c(55, 42)
# Total n obsevations for t = 1
n_tox1 <- c(165, 118)
# Cases for t = 1
r_tox2 <- c(77, 61, 61, 25)
# Total n obsevations for t = 2
n_tox2 <- c(161, 203, 143, 46)

# Data on time to response (sub-model 3) ----------------------------------
# Observed median time
x_tr <- c(23, 12)
# Observed variance
phi_tr <- c(3, 2) ^ 2

# Data on time to progression (sub-model 4) -------------------------------
# Number of studies
N_tp <- c(2, 3)
# Observed median time for t = 1
x_tp1 <- c(23, 21)
# Observed variance for t = 1
phi_tp1 <- c(4.077, 1.943) ^ 2
# Observed median time for t = 2
x_tp2 <- c(19, 26, 27.3)
# Observed variance for t = 2
phi_tp2 <- c(2.048, 2.263, 1.631) ^ 2
# Markov Model cycle length
tau <- 3

# Data on survival time (sub-model 5) -------------------------------------
# Number of studies
N_surv <- c(2, 3)
# Observed median time for t = 1
x_surv1 <- c(60.66, 47)
# Observed variance for t = 1
phi_surv1 <- c(3.9723, 5.6880) ^ 2
# Observed median time for t = 2
x_surv2 <- c(47.67, 65, 45.07)
# Observed variance for t = 2
phi_surv2 <- c(4.845, 5.171, 2.174) ^ 2
# Computed precision for t = 2

# JAGS model --------------------------------------------------------------
model_string <- "model {
 # Sub-model 1
 for (i in 1:N_res[1]) {
  r_res1[i] ~ dbin(p_res1[i], n_res1[i])
  logit(p_res1[i]) <- delta_res1[i]
  delta_res1[i] ~ dnorm(mu_res[1], tau_res[1])
 }
 for (i in 1:N_res[2]) {
  r_res2[i] ~ dbin(p_res2[i], n_res2[i])
  logit(p_res2[i]) <- delta_res2[i]
  delta_res2[i] ~ dnorm(mu_res[2], tau_res[2])
 }
 for (t in 1:2) {
  mu_res[t] ~ dnorm(0, 10)
  sigma_res[t] ~ dunif(0, 100)
  tau_res[t] <- 1 / (sigma_res[t] * sigma_res[t])
  pi_res[t] <- exp(mu_res[t]) / (1 + exp(mu_res[t]))
 }
  # Sub-model 2
 for (i in 1:N_tox[1]) {
  r_tox1[i] ~ dbin(p_tox[i], n_tox1[i])
  logit(p_tox[i]) <- delta_tox1[i]
  delta_tox1[i] ~ dnorm(mu_tox[1], tau_tox[1])
 }
 for (i in 1:N_tox[2]) {
  r_tox2[i] ~ dbin(p_tox2[i], n_tox2[i])
  logit(p_tox2[i]) <- delta_tox2[i]
  delta_tox2[i] ~ dnorm(mu_tox[2], tau_tox[2])
 }
 for (t in 1:2) {
  mu_tox[t] ~ dnorm(0, 10)
  sigma_tox[t] ~ dunif(0, 100)
  tau_tox[t] <- 1 / (sigma_tox[t] * sigma_tox[t])
  pi_tox[t] <- exp(mu_tox[t]) / (1 + exp(mu_tox[t]))
 }
  # Sub-model 3
 for (t in 1:2) {
  x_tr[t] ~ dnorm(mu_tr[t], prec_tr[t])
  mu_tr[t] ~ dnorm(0, 10)T(0, )
  rho_tr[t] <- -log(.5) / mu_tr[t]
  beta_tr[t] <- 1 - exp(-rho_tr[t] / tau)
  prec_tr[t] <- 1 / phi_tr[t]
 }
  # Sub-model 4
 for (i in 1:N_tp[1]) {
  x_tp1[i] ~ dnorm(delta_tp1[i], prec_tp1[i])
  delta_tp1[i] ~ dnorm(mu_tp[1], prec_mu_tp[1])
  prec_tp1[i] <- 1 / phi_tp1[i]
 }
 for (i in 1:N_tp[2]) {
  x_tp2[i] ~ dnorm(delta_tp2[i], prec_tp2[i])
  delta_tp2[i] ~ dnorm(mu_tp[2], prec_mu_tp[2])
  prec_tp2[i] <- 1 / phi_tp2[i]
 }
 for (t in 1:2) {
  mu_tp[t] ~ dnorm(0, prec_mu_tp[t])T(0, )
  sigma_tp[t] ~ dunif(0, 100)
  prec_mu_tp[t] <- 1 / (sigma_tp[t] * sigma_tp[t])
  rho_tp[t] <- -log(.5) / mu_tp[t]
  beta_tp[t] <- 1 - exp(-rho_tp[t] / tau)
 }
  # Sub-model 5
 for (i in 1:N_surv[1]) {
  x_surv1[i] ~ dnorm(delta_surv1[i], prec_surv1[i])
  delta_surv1[i] ~ dnorm(mu_surv[1], prec_mu_surv[1])
  prec_surv1[i] <- 1 / phi_surv1[i]
 }
 for (i in 1:N_surv[2]) {
  x_surv2[i] ~ dnorm(delta_surv2[i], prec_surv2[i])
  delta_surv2[i] ~ dnorm(mu_surv[2], prec_mu_surv[2])
  prec_surv2[i] <- 1 / phi_surv2[i]
 }
 for (t in 1:2) {
  mu_surv[t] ~ dnorm(0, prec_mu_surv[t])T(0, )
  sigma_surv[t] ~ dunif(0, 10)
  prec_mu_surv[t] <- 1 / (sigma_surv[t] * sigma_surv[t])
  rho_surv[t] <- -log(.5) / mu_surv[t]
  beta_surv[t] <- 1 - exp(-rho_surv[t] / tau)
 }
 for (t in 1:2) {
 mu_dth[t] <- mu_surv[t] - mu_tp[t]
 }
}
"
writeLines(text = model_string, con = "EvSynthMMExample.txt")

data_JAGS <- c("N_res", "n_res1", "n_res2", "N_surv", "N_tox", "n_tox1", "n_tox2",
               "N_tp", "phi_surv1", "phi_surv2", "phi_tp1", "phi_tp2", "phi_tr",
               "r_res1", "r_res2", "r_tox1", "r_tox2", "tau", "x_surv1", "x_surv2",
               "x_tp1", "x_tp2", "x_tr")
data_JAGS

params <- c("pi_res", "pi_tox", "beta_tr", "beta_tp", "beta_surv", "mu_dth")
jags_mod <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "EvSynthMMExample.txt", n.iter = 40000, n.chains = 4,
                 n.thin = 80)
jags_mod
attach.jags(jags_mod)

# Visual inspection of posterior ------------------------------------------
posterior <- as.array(jags_mod$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("pi_res[1]", "pi_res[2]",
                               "pi_tox[1]", "pi_tox[2]", "beta_tr[1]", "beta_tr[2]",
                               "beta_tp[1]", "beta_tp[2]", "beta_surv[1]", "beta_surv[2]"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , 1:13], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("pi_res[1]", "pi_res[2]",
                               "pi_tox[1]", "pi_tox[2]", "beta_tr[1]", "beta_tr[2]",
                               "beta_tp[1]", "beta_tp[2]", "beta_surv[1]", "beta_surv[2]"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("pi_res[1]", "pi_tox[1]", "beta_tr[1]", "beta_tp[1]", 
                               "beta_surv[1]"),
           off_diag_args = list(size = 1.5))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("pi_res[2]", "pi_tox[2]", "beta_tr[2]",
                               "beta_tp[2]", "beta_surv[2]"),
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("pi_res[1]", "pi_res[2]",
                               "pi_tox[1]", "pi_tox[2]", "beta_tr[1]", "beta_tr[2]",
                               "beta_tp[1]", "beta_tp[2]", "beta_surv[1]", "beta_surv[2]"), 
         lags = 150)

# Misc: define treatments ---------------------------------------------
colnames(pi_res) <- c("Status Quo", "New Treatment")
colnames(pi_res)
pi_res
colnames(pi_tox) <- c("Status Quo", "New Treatment")
colnames(pi_tox)
colnames(beta_tr) <- c("Status Quo", "New Treatment")
colnames(beta_tr)
colnames(beta_tp) <- c("Status Quo", "New Treatment")
colnames(beta_tp)
colnames(beta_surv) <- c("Status Quo", "New Treatment")
colnames(beta_surv)
colnames(mu_dth) <- c("Status Quo", "New Treatment")
colnames(mu_dth)

# Compute rho_dth and beta_dth values -------------------------------------
rho_dth <- matrix(NA, dim(mu_dth), ncol = 2)

rho_dth[, 1] <- -log(.5) / mu_dth[, 1]
rho_dth[, 2] <- -log(.5) / mu_dth[, 2]
rho_dth

beta_dth <- matrix(NA, dim(mu_dth), ncol = 2)
beta_dth[, 1] <- 1 - exp(-rho_dth[, 1] / 3)
beta_dth[, 2] <- 1 - exp(-rho_dth[, 2] / 3)
colnames(beta_dth) <- c("Status Quo", "New Treatment")
colnames(beta_dth)

# ==========================================================================================
# Markov Model setup: -----------------------------------------------------------
# ==========================================================================================
n_age_init <- 15 # age at baseline
n_age_max <- 100 # maximum age of follow up
n_t <- n_age_max - n_age_init # time horizon, number of cycles

v_n <- c("Stable", "Response", "Progression", "Death") # the 4 health states of the model:
n_states <- length(v_n) # number of health states 

d_c <- 0.03 # discount rate for costs 
d_e <- 0.03 # discount rate for QALYs
v_names_str <- c("Usual care", "New treatment") # store the strategy names

a_P <- array(0, dim = c(n_states, n_states, n_t, n.sims),
             dimnames = list(v_n, v_n, 0:(n_t - 1), 1:n.sims))
dim(a_P)
# NOTE: WHEN READY ADD TREATMENTS 1 AND 2, I.E. MAKE a_P_1 AND a_P_2... etc.

# ==========================================================================================
# Fill in transition array ------------------------------------------------
# ==========================================================================================
# Status Quo  ----------------------------------------------------------


# New Treatment  ----------------------------------------------------------
# All transitions from the state Stable
a_P["Stable", "Stable", , ] <- (1 - pi_tox[, 2]) * (1 - beta_tr[, 2])
a_P["Stable", "Progression", , ] <- pi_tox[, 2]
a_P["Stable", "Response", , ] <- (1 - pi_tox[, 2]) * beta_tr[, 2]
a_P["Stable", "Death", , ] <- 1 - ((1 - pi_tox[, 2]) * (1 - beta_tr[, 2]) + pi_tox[, 2] + 
                                    (1 - pi_tox[, 2]) * beta_tr[, 2])
# All transitions from the state Response
a_P["Response", "Response", , ] <- ((1 - pi_tox[, 2] ) * pi_res[, 2])
a_P["Response", "Progression", , ] <-  (pi_tox[, 2]) + (1 - pi_tox[, 2]) * (1 - pi_res[, 2])

# All transitions from the state Progression
a_P["Progression", "Progression", , ] <- 1 - beta_dth[, 2]
a_P["Progression", "Death", , ] <- beta_dth[, 2]

# All transitions from the state Death
a_P["Death", "Death", , ] <- 1
a_P

# ==========================================================================================
# Probability Validation ----------------------------------------------
# ==========================================================================================
# For each state, all transitions must add up to the total number of n_sims, or 1 when 
# divided by the number of n_sims.

# Stable:
sum((1 - pi_tox[, 2]) * (1 - beta_tr[, 2]) + (pi_tox[, 2]) + 
      ((1 - pi_tox[, 2]) * beta_tr[, 2]) + 1 - ((1 - pi_tox[, 2]) * (1 - beta_tr[, 2]) + 
                                                  pi_tox[, 2] + 
                                    (1 - pi_tox[, 2]) * beta_tr[, 2])) / n.sims
# Response:
sum(((1 - pi_tox[, 2] )* pi_res[, 2]) + 
      (pi_tox[, 2]) + (1 - pi_tox[, 2]) * (1 - pi_res[, 2])) / n.sims
# Progression:
sum((1 - beta_dth[, 2]) + beta_dth[, 2]) / n.sims
# and death is obvious...

# ==========================================================================================
# Run Markov model --------------------------------------------------------
# ==========================================================================================
## Initial state vector
# All starting healthy
v_s_init <- c(Stable = 1, Response = 0, Progression = 0, Death = 0) # initial state vector
v_s_init

## Initialize cohort trace for age-dependent cSTM
m_M_ad <- array(matrix(0, nrow = n_t + 1, ncol = n_states),
                dim = c(c(n_t + 1, n_states), n.sims), 
                dimnames = list(0:n_t, v_n, 1:n.sims))
dim(m_M_ad)
m_M_ad

# Store the initial state vector in the first row of the cohort trace
m_M_ad[1, , ] <- v_s_init
# Check:
m_M_ad[1:2, ,  ]

## Initialize transition array for each cycle, for each n.sim:
a_A <- array(0,
             dim = c(n_states, n_states, n_t + 1, n.sims),
             dimnames = list(v_n, v_n, 0:n_t, 1:n.sims))
dim(a_A)
# Set first slice of A with the initial state vector in its diagonal
for (i in 1:n.sims) {
 diag(a_A[, , 1, i]) <- v_s_init
}
# Check (first state, all states, for all cycles, across all n.sims):
a_A[1, 1:4, , ]

# Iterative solution of age-dependent cSTM
for (t in 1:n_t) {
 for(i in 1:n.sims) {
  # Fill in cohort trace
  m_M_ad[t + 1, , i] <- m_M_ad[t, , i] %*% a_P[ , , t, i]
  # Fill in transition dynamics array
  a_A[, , t + 1, i]  <- m_M_ad[t, , i] * a_P[ , , t, i]
  }
}
str(a_P)
str(m_M_ad)

# Visualise Markov Model -----------------------------------------------------
# code for the DARTH colors for the figures
DARTHgreen      <- '#9A031E'  
DARTHyellow     <- '#E36414'  
DARTHblue       <- '#FFD959' 
DARTHlightgreen <- '#D9CFC1'
DARTHgray       <- '#0F4C5C'

cols <- c("Stable" = DARTHgreen, "Response" = DARTHblue, 
          "Progression" = DARTHyellow, "Death" = DARTHgray)
lty <-  c("Stable" = 1, "Response" = 2, "Progression" = 3, "Death" = 4)

ggplot(melt(apply(m_M_ad, c(1, 2), mean)), aes(x = Var1, y = value, 
                      color = Var2, linetype = Var2)) +
 geom_line(size = 1) +
 scale_colour_manual(name = "Health state", 
                     values = cols) +
  scale_linetype_manual(name = "Health state",
                       values = lty) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  theme_light(base_size = 14) +
  theme(legend.position = "bottom", 
        legend.background = element_rect(fill = NA),
        text = element_text(size = 15))

barplot(apply(m_M_ad, c(2, 1), mean), space = 1,
        ylab = "Proportion of patients in each state",
        main = "New Treatment", col = cols, cex.names = 1, cex.main = 1.5)
legend("topright", legend = c("Stable", "Response", "Progrssion", "Death"), 
       fill = cols, cex = .5, box.lwd = 2,
       horiz = FALSE)

# Survival curve:
v_S_ad <- rowSums(m_M_ad[, -4, ]) / n.sims # vector with survival curve
v_S_ad
ggplot(data.frame(Cycle = 0:n_t, Survival = v_S_ad), 
       aes(x = Cycle, y = Survival)) +
  geom_line(size = 1.3) +
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = 14) +
  theme()

# Life expectancy
le_ad <- sum(v_S_ad)
le_ad

# Cost-effectiveness analysis ---------------------------------------------


