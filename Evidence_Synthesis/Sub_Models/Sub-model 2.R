library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(parallel)

options(mc.cores = detectCores())

# The following sub-model simulates age-specific probability of infection using informative 
# priors and collated data from several observational studies.

# ==========================================================================================
# Age-specific HPV incidence ----------------------------------------------------
# ==========================================================================================

# Informative prior -------------------------------------
# Study: Sinanovic, E., et al. 2009
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Age-specific incidence:
age_group <- c("15-16", "17", "18", "19", "20", "21", "22-23", "24-29", "30-49", "50≥")
gr_Age <- c(rep(age_group[1], length(15:16)), rep(age_group[2], 1), 
           rep(age_group[3], 1), rep(age_group[4], 1), rep(age_group[5], 1), 
           rep(age_group[6], 1), rep(age_group[7], length(22:23)), 
           rep(age_group[8], length(24:29)), rep(age_group[9], length(30:49)), 
           rep(age_group[10], length(50:85)))
gr_Age

incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, .1, .05, .01, .005))
rate_Age <- c(rep(incidence[1], length(15:16)), rep(incidence[2], 1), 
           rep(incidence[3], 1), rep(incidence[4], 1), rep(incidence[5], 1), 
           rep(incidence[6], 1), rep(incidence[7], length(22:23)), 
           rep(incidence[8], length(24:29)), rep(incidence[9], length(30:49)), 
           rep(incidence[10], length(50:85)))
rate_Age
# Convert to probability:
pPrior_Age <- 1 - exp(- rate_Age * 1)
pPrior_Age

barplot(rate_Age, names.arg = "From 15 to 90 years of age", 
        main = "Rate of HPV+", sub =  "At constant t-cycle of 1-year")

# Study: Peto J., et al 2004 --------------------------------------------
# Cervical HPV infection and neoplasia in a large population-based prospective study: the 
# Manchester cohort.

# 15-19:
# n = 319
# r = 69

# 20-24:
# n = 466
# r = 92

# 25-29:
# n = 586
# r = 93

# 30-34:
# n = 1191
# r = 86

# 35-39:
# n = 914
# r = 47

# 40-44:
# n = 926
# r = 28

# 45-49:
# n = 508
# r = 15

# 50-54:
# n = 1104
# r = 29

# 55-59:
# n = 114
# r = 1

n_Age <- c(rep(319, length(15:19)), rep(466, length(20:24)), rep(586, length(25:29)), 
           rep(1191, length(30:34)), rep(914, length(35:39)), rep(926, length(40:44)), 
           rep(508, length(45:49)), rep(1104, length(50:54)), rep(114, length(55:85))
           )
n_Age

r_Age <- c(rep(69, length(15:19)), rep(92, length(20:24)), rep(93, length(25:29)), 
           rep(86, length(30:34)), rep(47, length(35:39)), rep(28, length(40:44)), 
           rep(15, length(45:49)), rep(29, length(50:54)), rep(1, length(55:85))
           )
r_Age

# MM time horizon:
n_t <- length(r_Age)
n_t

# Model -------------------------------------------------------------------
model_String <- "
 model {
  for (i in 1:n_t) {
   # Likelihood:
   r_Age[i] ~ dbin(omega_Age[i], n_Age[i])
   # Logit function:
   logit(omega_Age[i]) <- nu[i] + eta[i]
  
  # Random effect:
  eta[i] ~ dnorm(0, 1.0e-5)
  # Informative prior:
  nu[i] ~ dnorm(iota, phi.prec)
  
  }
  # Hyperpriors on precision:
  iota ~ dnorm(0, 1.0e-6)
  phi ~ dunif(0, 10)
  phi.prec <- 1 / (phi * phi)
  
}
"
writeLines(text = model_String, con = "Inc_Prob.txt")

data_JAGS <- list(
 r_Age = r_Age,
 n_Age = n_Age,
 n_t = n_t
)
data_JAGS

inits <- list(
 list(phi = .1, eta = c(rep(0, n_t))),
 list(phi = .1, eta = c(rep(0, n_t)))
 )

params <- c("omega_Age")

mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "Inc_Prob.txt", inits = inits,
                 n.chains = 2, n.iter = 10000, n.burnin = 1000, 
                 n.thin = 20)
mod_JAGS

# Attach JAGS model to global envir:
attach.jags(mod_JAGS)

# Create array of posterior samples:
posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

# Inspect trace of chain for each parameter:
color_scheme_set("viridisD")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("omega_Age[1]", "omega_Age[4]", "omega_Age[6]", 
                               "omega_Age[70]", "deviance"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , 1:10], window = c(341, 431), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("omega_Age[1]", "omega_Age[4]", "omega_Age[6]", 
                               "omega_Age[9]", "deviance"))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("omega_Age[1]", "omega_Age[4]", "omega_Age[6]", 
                               "omega_Age[9]", "deviance"),
         lags = 250)

# Calibration -------------------------------------------------------------
barplot(rate_Age, axes = FALSE)
par(new = TRUE)
plot(apply(omega_Age, 2, mean), type = "l", col = "red", lwd = 2,
     ylab = "Probability of infection", xlab = "By Age: 15 - 85 years", x = 15:85, 
     xgap.axis = 15)

# End file ----------------------------------------------------------------