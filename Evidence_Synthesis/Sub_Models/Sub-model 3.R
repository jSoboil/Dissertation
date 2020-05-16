library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(parallel)

options(mc.cores = detectCores())

# The following sub-model simulates the average vaccine efficiacy using data from multiple 
# RCTs.

# ===========================================================================================
# Vaccine efficacy --------------------------------------------------------
# ===========================================================================================
# Vaccine efficacy is expressed as the reduction in the occurrence of HPV due to vaccination.
# Individuals who are vaccinated fully experience a rate of occurrence of HPV that is 
# (1 – alpha) times that of those who are not vaccinated (Favato G., et al. 2012).

# Study 1: Munoz N., et al. 2009 ------------------------------------------
# Safety, immunogenicity, and eðcacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) recombinant vaccine in women aged 24–45 years: a randomised, double-blind trial.

# Women aged 24–45 years with no history of genital warts or cervical disease were enrolled 
# from community health centres, academic health centres, and primary health-care providers 
# into an ongoing multicentre, parallel, randomised, placebo-controlled, double-blind study.

# Events by group:
# tB:4
# nB:1615

# tA:41
# nA:1607

 # Study 2: Villa L, et al. 2006 -------------------------------------------
# High sustained efficacy of a prophylactic quadrivalent human papillomavirus types 
# 6/11/16/18 L1 virus-like particle vaccine through 5 years of follow-up.

# The study enrolled nonpregnant, healthy women who had no prior abnormal Pap smears, and 
# reported a lifetime history of four or fewer male sex partners.

# Events by group:
# tB:2
# nB:235

# tA:45
# nA:233

# Study 3: Catellsague X., et al. 2011 ------------------------------------
# End-of-study safety, immunogenicity, and efficacy of quadrivalent HPV (types 6, 11, 16, 
# 18) recombinant vaccine in adult women 24–45 years of age.

# Study enrolled 3819 24–45-year-old women with no history of cervical disease or genital 
# warts in the past 5 years.

# Events by group:
# tB:1
# nB:1578

# tA:38
# nA:1583

# Study 4: Wei L., et al. 2018 --------------------------------------------
# Efficacy of quadrivalent human papillomavirus vaccine against persistent infection and 
# genital disease in Chinese women: A randomized, placebo-controlled trial with 78-month 
# follow-up.

# Pregnant women and those with a history of genital warts or significant cervical disease, 
# active cervical disease, or prior HPV vaccine recipients were excluded.

# Events by group:
# tB:3
# nB:1271

# tA:48
# nA:1243


# Study 5: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tB:3
# nB:1990

# tA:25
# nA:1880

# Study 6: Konno R., et al. 2010 ------------------------------------------
# Efficacy of Human Papillomavirus Type 16/18 AS04YAdjuvanted Vaccine in Japanese Women Aged
# 20 to 25 Years.

# Final Analysis of a Phase 2 Double-Blind, Randomized Controlled Trial

# Events by group:
# tB:0
# nB:387

# tA:15
# nA:392

# Study 7: Yoshikawa H., et al. 2013 --------------------------------------
# Efficacy of quadrivalent human papillomavirus (types 6, 11, 16 and 18) vaccine (GARDASIL) 
# in Japanese women aged 18–26 years.

# A randomized double-blind placebo-controlled phase II trial was conducted to evaluate the 
# efficacy of a prophylactic quadrivalent vaccine targeting the human papillomavirus (HPV) 
# types most fre- quently associated with cervical cancer (types 16 ⁄ 18) and genital warts 
# (types 6 ⁄ 11) in Japanese women aged 18–26 years.

# Events by group:
# tB:3
# nB:419

# tA:24
# nA:422

# Vaccine efficacy data ----------------------------------------------------
# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 7
rA.vac <- c(41, 45, 38, 48, 25, 15, 24)
nA.vac <- c(1607, 233, 1583, 1243, 1880, 392, 422)
rB.vac <- c(4, 2, 1, 3, 3, 0, 3)
nB.vac <- c(1615, 235, 1578, 1271, 1990, 387, 419)

# ==========================================================================================
# Sub-models 3 ------------------------------------------------------
# ==========================================================================================
model_String <- "model {
# Sub-model 3: vaccine-efficacy against infection:
  # model parmaeters abbreviated by .vac
 for (i in 1:Nstud.vac) {
  rA.vac[i] ~ dbin(pA.vac[i], nA.vac[i])
  rB.vac[i] ~ dbin(pB.vac[i], nB.vac[i])
 
  logit(pA.vac[i]) <- mu.vac[i]
  logit(pB.vac[i]) <- mu.vac[i] + delta.vac[i]
 
  mu.vac[i] ~ dnorm(0, 1.0e-5)
  delta.vac[i] ~ dnorm(psi.vac, prec.vac)
}
# Priors for sub-model 3:
 psi.vac ~ dnorm(0, 1.0e-6)
 tau.vac ~ dunif(0, 10)
 prec.vac <- 1 / pow(tau.vac, 2)
 
 OR.vac <- exp(psi.vac)
 pEfficacy.vac <- (OR.vac / 1 + OR.vac)
 
}
"
writeLines(text = model_String, con = "Vac_Efficacy.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(Nstud.vac = Nstud.vac, rA.vac = rA.vac, rB.vac = rB.vac, nA.vac = nA.vac, 
                  nB.vac = nB.vac)
# Initial JAGS sampler values:
# Parameters to monitor:
params <- c("OR.vac", "pEfficacy.vac")

mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "Vac_Efficacy.txt", n.chains = 2, 
                 n.iter = 10000, n.burnin = 1000, n.thin = 18)
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
         lags = 250)

# Vaccine efficacy probability:
 1 - apply(pEfficacy.vac, 2, mean)
# Example: say 21% of an age pop. has HPV:
.21 * apply(pEfficacy.vac, 2, mean)

# End file ----------------------------------------------------------------