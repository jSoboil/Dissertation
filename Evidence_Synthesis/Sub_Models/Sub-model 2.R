library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(parallel)

options(mc.cores = detectCores())

# The following sub-model simulates age-specific prevalence and average vaccine efficiacy 
# using data from multiple RCTs. Assumed any end-point - infection or disease - signalling 
# efficacy.

# ==========================================================================================
# Age-specific infection ----------------------------------------------------
# ==========================================================================================
# Note: will use a relative risk parameter for proportion of HIV+ population who have 
# increased risk of being reinfected rather than clearing HPV. See Italian study.

# Informative prior -------------------------------------
# Data collected from HPV and Related Diseases Report:
age_group <- c("≤25", "25-34", "35-44", "45-54", "55-64")
# Note: probability converted to yearly rate for use in lognormal distribution.
# r = - (1 / t) * log(1 - p)
incidence <- c(0.5798185, 0.4620355, 0.2231436, 0.1984509, 0.1863296)

mu.log <- lnorm_params(m = incidence, v = 1)$mu
sigma.log <- lnorm_params(m = incidence, v = 1)$sigma
mu.log
sigma.log
cbind(age_group, incidence)

# Prior sampling Model :
# omega.age[i] ~ dlnorm(mu.log[i], sigma.log[i])

# ===========================================================================================
# Vaccine efficacy --------------------------------------------------------
# ===========================================================================================
# Vaccine efficacy is expressed as the reduction in the occurrence of HPV due to vaccination.
# Individuals who are vaccinated fully experience a rate of occurrence of HPV that is 
# (1 – alpha) times that of those who are not vaccinated (Favato G., et al. 2012).

# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B).

# Study 1: Munoz N., et al. 2009 ------------------------------------------
# Safety, immunogenicity, and eðcacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) recombinant vaccine in women aged 24–45 years: a randomised, double-blind trial.

# Women aged 24–45 years with no history of genital warts or cervical disease were enrolled 
# from community health centres, academic health centres, and primary health-care providers 
# into an ongoing multicentre, parallel, randomised, placebo-controlled, double-blind study.

# Efficacy against the incidence of infection .

# Events by group:
# tA:41
# nA:1607
# tB:4
# nB:1615

 # Study 2: Villa L, et al. 2006 -------------------------------------------
# High sustained efficacy of a prophylactic quadrivalent human papillomavirus types 
# 6/11/16/18 L1 virus-like particle vaccine through 5 years of follow-up.

# The study enrolled nonpregnant, healthy women who had no prior abnormal Pap smears, and 
# reported a lifetime history of four or fewer male sex partners. By end point persistent 
# infection.

# Events by group:
# tA:46
# nA:233
# tB:2
# nB:235

# Study 3: Catellsague X., et al. 2011 ------------------------------------
# End-of-study safety, immunogenicity, and efficacy of quadrivalent HPV (types 6, 11, 16, 
# 18) recombinant vaccine in adult women 24–45 years of age.

# Study enrolled 3819 24–45-year-old women with no history of cervical disease or genital 
# warts in the past 5 years. By end point persistent infection.

# Events by group:
# tA:38
# nA:1583
# tB:1
# nB:1578

# Study 4: Wei L., et al. 2018 --------------------------------------------
# Efficacy of quadrivalent human papillomavirus vaccine against persistent infection and 
# genital disease in Chinese women: A randomized, placebo-controlled trial with 78-month 
# follow-up.

# Pregnant women and those with a history of genital warts or significant cervical disease, 
# active cervical disease, or prior HPV vaccine recipients were excluded. 12-Month cervical
# Persisitent Infection.

# Events by group:
# tA:53
# nA:1245
# tB:5
# nB:1276

# Study 5: Yoshikawa H., et al. 2013 --------------------------------------
# Efficacy of quadrivalent human papillomavirus (types 6, 11, 16 and 18) vaccine (GARDASIL) 
# in Japanese women aged 18–26 years.

# A randomized double-blind placebo-controlled phase II trial was conducted to evaluate the 
# efficacy of a prophylactic quadrivalent vaccine targeting the human papillomavirus (HPV) 
# types most fre- quently associated with cervical cancer (types 16 ⁄ 18) and genital warts 
# (types 6 ⁄ 11) in Japanese women aged 18–26 years.

# Events by group:
# tA:24
# nA:422
# tB:3
# nB:419

# Study 6: Garland S., et al. Future I. 2007 --------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# A phase 3 trial was conducted to evaluate the efficacy of a prophylactic quadrivalent 
# vaccine in preventing anogenital diseases associated with human papillomavirus (HPV) types
# 6, 11, 16, and 18.

# Events by group:
# tA:60
# nA:2279
# tB:0
# nB:2261

# Study 7: Koutsky L., et al. Future II. 2007 -----------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent High-Grade Cervical Lesions.

# Randomized, double-blind trial, we assigned 12,167 women between the ages of 15 and 26 
# years to receive three doses of either HPV-6/11/16/18 vaccine or placebo, administered at 
# day 1, month 2, and month 6. The primary analysis was performed for a per-protocol 
# susceptible population that included 5305 women in the vaccine group and 5260 in the 
# placebo group who had no virologic evidence of infection with HPV-16 or HPV-18 through 1 
# month after the third dose (month 7).

# Events by group:
# tA:42
# nA:5260
# tB:1
# nB:5305

# Vaccine efficacy data ----------------------------------------------------
# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 7
rA.vac <- c(41, 46, 38, 53, 24, 60, 42)
nA.vac <- c(1607, 233, 1583, 1245, 422, 2279, 5260)
rB.vac <- c(4, 2, 1, 5, 3, 0, 1)
nB.vac <- c(1615, 235, 1578, 1276, 419, 2261, 5305)

# ==========================================================================================
# Sub-models 3 ------------------------------------------------------
# ==========================================================================================
model_String <- "model {
  # Sub-model 2: age-specific prevalence:
  for (i in 1:5) {
    omega.age[i] ~ dlnorm(mu.log[i], sigma.log[i])
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
  psi.vac ~ dnorm(0, 1)
  tau.vac ~ dunif(0, 1)
  prec.vac <- 1 / pow(tau.vac, 2)
  
  OR.vac <- exp(psi.vac)
  pEfficacy.vac <- (OR.vac / 1 + OR.vac)
 
}
"
writeLines(text = model_String, con = "Age_and_Efficacy.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(Nstud.vac = Nstud.vac, rA.vac = rA.vac, rB.vac = rB.vac, 
                  nA.vac = nA.vac, nB.vac = nB.vac, mu.log = mu.log, 
                  sigma.log = sigma.log)
# Initial JAGS sampler values:
# Parameters to monitor:
params <- c("OR.vac", "pEfficacy.vac", "omega.age")
# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 20000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 250)

mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "Age_and_Efficacy.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS
# Attach JAGS model to global envir:
attach.jags(mod_JAGS)

# Create array of posterior samples:
posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

# Inspect trace of chain for each parameter:
color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("pEfficacy.vac", "omega.age[1]", "omega.age[2]", 
                               "omega.age[3]", "omega.age[4]", "omega.age[5]"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , ], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("pEfficacy.vac", "omega.age[1]", "omega.age[2]", 
                               "omega.age[3]", "omega.age[4]", "omega.age[5]"))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("pEfficacy.vac", "omega.age[1]", "omega.age[2]", 
                               "omega.age[3]", "omega.age[4]", "omega.age[5]"),
         lags = 150)

# Calibration: ------------------------------------------------------------
# Age-specific prevalence:
plot(apply(omega.age, 2, mean), type = "l", ylab = "Prevalence", xlab = "By Age", 
     xgap.axis = 100)
# Vaccine efficacy probability:
1 - apply(pEfficacy.vac, 2, mean)
# Example: say 21% of some age pop. group has HPV:
.21 * apply(pEfficacy.vac, 2, mean)
# End file ----------------------------------------------------------------