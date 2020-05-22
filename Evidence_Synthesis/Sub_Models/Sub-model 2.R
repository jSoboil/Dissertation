library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(parallel)

options(mc.cores = detectCores())

# The following sub-model simulates the risk increase that faces the HIV+ proportion of 
# individuals in South Africa.

# ===========================================================================================
# Relative Risk Increase (HIV+) ----------------------------------------------------
# ===========================================================================================

# Study: Mbulawa Z., et al. 2015 ------------------------------------------
# Human papillomavirus prevalence in South African women and men according to age and human 
# immunodeficiency virus status.

# Then use proportion of HIV+ pop. for how many at risk at time j. See Italian study.

# Participants were 208 HIV-negative women, 278 HIV-positive women, 325 HIV-negative men and
# 161 HIV-positive men between the ages of 18–66 years. HPV types were determined in cervical
# and penile cells by Roche Linear Array HPV genotyping assay.

# HIV+ women positive to HPV: 205/277
# HIV- women positive to HPV: 76/207

rHIV_pos <- 205
nHIV_pos <- 277
rHIV_neg <- 76
nHIV_neg <- 207

# Proportion of pop. at increased risk due to HIV:
zeta_alpha <- beta_params(mean = .131, sigma = .05)$alpha
zeta_beta <- beta_params(mean = .131, sigma = .05)$beta

model_String <- "model {
# SUB-MODEL 4: POPULATION AT INCREASED RISK 
# DUE TO HIV

# Binomial Likelihood:
 rHIV_pos ~ dbin(pHIV_pos, nHIV_pos)
 rHIV_neg ~ dbin(pHIV_neg, nHIV_neg)
 # Logistic function:
 logit(pHIV_pos) <- mu.hiv + delta.hiv
 logit(pHIV_neg) <- mu.hiv
 
 # Average effect:
 mu.hiv ~ dnorm(0, 1 / 10)
 # Random effect:
 delta.hiv ~ dnorm(nu.hiv, prec.hiv)
 
 # Hyperpriors on delta.hiv:
 nu.hiv ~ dnorm(0, 1 / 10)
 prec.hiv <- 1 / psi.hiv * psi.hiv
 psi.hiv ~ dunif(0, 1)
 
 HIV_odds <- exp(nu.hiv)
 HIV_risk <- psi.hiv / (1 + psi.hiv)
 
 # Proportion of population at
 # increased risk:
 zeta_HIVpos ~ dbeta(zeta_alpha, zeta_beta)
 
}
"
writeLines(text = model_String, con = "HIV_risk.txt")

data_JAGS <- list(rHIV_pos = rHIV_pos, rHIV_neg = rHIV_neg, 
                  nHIV_pos = nHIV_pos, nHIV_neg = nHIV_neg,
                  zeta_alpha = zeta_alpha, zeta_beta = zeta_beta)

params <- c("HIV_odds", "HIV_risk", "zeta_HIVpos")

mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "HIV_risk.txt", n.iter = 20000, n.burnin = 5000, 
                 n.thin = 44)
mod_JAGS
# Attach JAGS model to global envir:
attach.jags(mod_JAGS)

# Create array of posterior samples:
posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

# Inspect trace of chain for each parameter:
color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = params,
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , ], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = params)

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = params,
         lags = 150)

# End file ----------------------------------------------------------------