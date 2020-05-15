library(tidyverse)
library(rjags)
library(R2jags)
library(bayesplot)
library(parallel)

options(mc.cores = detectCores())

# The following sub-models simulate age-specific probability of infection using informative 
# priors, and the average vaccine efficiacy using data from multiple RCTs.

# ==========================================================================================
# Age-specific HPV incidence ----------------------------------------------------
# ==========================================================================================
# Study 1: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Age-specific incidence:
age_group <- c("15-16", "17", "18", "19", "20", "21", "22-23", "24-29", "30-49", "50≥")
group_Age <- c(rep(age_group[1], length(15:16)), rep(age_group[2], 1), 
           rep(age_group[3], 1), rep(age_group[4], 1), rep(age_group[5], 1), 
           rep(age_group[6], 1), rep(age_group[7], length(22:23)), 
           rep(age_group[8], length(24:29)), rep(age_group[9], length(30:49)), 
           rep(age_group[10], length(50:100)))
group_Age

incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, .1, .05, .01, .005))
r_Age <- c(rep(incidence[1], length(15:16)), rep(incidence[2], 1), 
           rep(incidence[3], 1), rep(incidence[4], 1), rep(incidence[5], 1), 
           rep(incidence[6], 1), rep(incidence[7], length(22:23)), 
           rep(incidence[8], length(24:29)), rep(incidence[9], length(30:49)), 
           rep(incidence[10], length(50:100)))
r_Age
barplot(r_Age, names.arg = "From Age 15 to 100", 
        main = "Rate of HPV+", sub =  "At constant t-cycle of 1-year")

# Function to convert means and st.deviations to the log-normal scale. Note: Copyright 
# Baio G. 2012.
lognPar <- function(mu, sdev) {
	sdev.sq <- sdev ^ 2
	mu.log <- log(mu) - .5 * log(1 + sdev.sq / mu ^ 2)
	sdev.sq.log <- log(1 + (sdev.sq / mu ^ 2))
	sigma.log <- sqrt(sdev.sq.log)
	list(mu.log = mu.log, sigma.log= sigma.log)
}

# Used as age-specific means for informative log-normal prior:
r_mu.log <- lognPar(mu = r_Age, sdev = r_Age)$mu.log
r_mu.log
# Used as age-specific st.deviation for informative log-normal prior:
r_sigma.log <- 1 / lognPar(mu = r_Age, sdev = r_Age)$sigma.log ^ 2
r_sigma.log
# Note: above computes the precision for the log-normal distribution outside of JAGS to speed
# up computation.

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

# Study 4: Garland S., et al. FUTURE I 2007 -------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# PPE pop:
# Events by group:
# tB:0
# nB:2261

# tA:60
# nA:2279

# Study 5: Wei L., et al. 2018 --------------------------------------------
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

# Study 6: Koutsky L., et al. FUTURE II Group. 2007 -----------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent High-Grade Cervical Lesions.

# PPE pop:
# Events by group:
# tB:1
# nB:5305

# tA:42
# nA:5260

# Study 7: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tB:3
# nB:1990

# tA:25
# nA:1880

# Vaccine efficacy data ----------------------------------------------------
# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 7
rA.vac <- c(41, 45, 38, 60, 48, 42, 25)
nA.vac <- c(1607, 233, 1583, 2279, 1243, 5260, 1880)
rB.vac <- c(4, 2, 1, 0, 3, 1, 3)
nB.vac <- c(1615, 235, 1578, 2261, 1271, 5305, 1990)

# ==========================================================================================
# Sub-models 2 and 3 ------------------------------------------------------
# ==========================================================================================
model_String <- "model {
# Sub-model 2: age-specific probability of infection:
  # model parameters abbreviated by .age
    
    # Note: equivalent of sampling directly from an
    # informative prior.
  for (i in 1:86) {
    omega.age[i] ~ dlnorm(r_mu.log[i], r_sigma.log[i])
  }

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
 tau.sq.vac <- (tau.vac * tau.vac)
 prec.vac <- 1 / tau.sq.vac
 
 OR.vac <- exp(psi.vac)
 pEfficacy.vac <- (OR.vac / 1 + OR.vac)
 
}
"
writeLines(text = model_String, con = "Inc_and_efficacy.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(Nstud.vac = Nstud.vac, rA.vac = rA.vac, rB.vac = rB.vac, nA.vac = nA.vac, 
                  nB.vac = nB.vac, r_mu.log = r_mu.log, 
                  r_sigma.log = r_sigma.log)
# Initial JAGS sampler values:
inits <- list(
 list(psi.vac = 1, tau.vac = .3, delta.vac = c(0, 0, 0, 0, 0, 0, 0), 
      mu.vac = c(0, 0, 0, 0, 0, 0, 0), 
      omega.age = c(rep(0, 30), rep(1, 30), rep(.5, 26))),
 list(psi.vac = 0, tau.vac = 1, delta.vac = c(0, 0, 0, 0, 0, 0, 0), 
      mu.vac = c(0, 0, 0, 0, 0, 0, 0),
      omega.age = c(rep(0, 30), rep(1, 30), rep(.5, 26)))
)
# Parameters to monitor:
params <- c("OR.vac", "pEfficacy.vac", "omega.age")

mod_JAGS <- jags(data = data_JAGS, inits = inits, parameters.to.save = params, 
                 model.file = "Inc_and_efficacy.txt", n.chains = 2, 
                 n.iter = 20000, n.burnin = 10000, n.thin = 20)
mod_JAGS
# Attach JAGS model to global envir:
attach.jags(mod_JAGS)

# Create array of posterior samples:
posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

# Inspect trace of chain for each parameter:
color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("OR.vac", "omega.age[1]", "omega.age[5]"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , 80:86], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("OR.vac", "omega.age[1]", "omega.age[5]"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("OR.vac", "omega.age[1]", "omega.age[5]"),
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("OR.vac", "omega.age[1]", "omega.age[5]"), 
         lags = 50)

# Convert rate of infection to probability of infection for Status Quo cohort:
p_Age <- 1 - exp(-r_Age * 1)
p_Age

# Vaccine efficacy probability:
1 - apply(pEfficacy.vac, 2, mean)
# Status quo compared to age dependet probability of infection given vaccine:
cbind("Status Quo" = p_Age, "Vaccine" = p_Age * apply(pEfficacy.vac, 2, mean))

# End file ----------------------------------------------------------------