library(rjags)
library(R2jags)
library(bayesplot)
library(parallel)

setwd(dir = "/Users/joshuamusson/Desktop/Analytics/R/Intergrated-CEA-thesis/Evidence_Synthesis/Sub_Models/")
getwd()
options(mc.cores = detectCores())

# ===========================================================================================
# Vaccine efficacy --------------------------------------------------------
# ===========================================================================================
# Vaccine efficacy, expressed as the reduction in the occurrence of HPV due to vaccination.
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

# ===========================================================================================
# Model specifications ----------------------------------------------------
# ===========================================================================================
# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B).

Nstud <- 7
rA <- c(41, 45, 38, 60, 48, 42, 25)
nA <- c(1607, 233, 1583, 2279, 1243, 5260, 1880)
rB <- c(4, 2, 1, 0, 3, 1, 3)
nB <- c(1615, 235, 1578, 2261, 1271, 5305, 1990)

model_String <- "model {
 for (i in 1:Nstud) {
  rA[i] ~ dbin(pA[i], nA[i])
  rB[i] ~ dbin(pB[i], nB[i])
 
  logit(pA[i]) <- mu[i]
  logit(pB[i]) <- mu[i] + delta[i]
 
  mu[i] ~ dnorm(0, 1.0e-5)
  delta[i] ~ dnorm(d, prec)
  
 }
 # Priors
 d ~ dnorm(0, 1.0e-6)
 tau ~ dunif(0, 10)
 tau.sq <- (tau * tau)
 prec <- 1 / tau.sq
 
 OR <- exp(d)
 pEfficacy <- OR / 1 + OR
}
"
writeLines(text = model_String, con = "VaccineEfficacy.txt")

# Transform data into format that can be read by JAGS:
data_JAGS <- list(Nstud = Nstud, rA = rA, rB = rB, nA = nA, nB = nB)
# Initial sampler values:
inits <- list(
 list(d = 0, tau = 1, delta = c(0, 0, 0, 0, 0, 0, 0), 
      mu = c(0, 0, 0, 0, 0, 0, 0)),
 list(d = 0, tau = 1, delta = c(0, 0, 0, 0, 0, 0, 0), 
      mu = c(0, 0, 0, 0, 0, 0, 0))
)
# Parameters to monitor:
params <- c("d", "OR", "pEfficacy")

mod_JAGS <- jags(data = data_JAGS, inits = inits, parameters.to.save = params, 
                 model.file = "VaccineEfficacy.txt", n.chains = 2, 
                 n.iter = 10000, n.burnin = 5000, n.thin = 10)
mod_JAGS
attach.jags(mod_JAGS)

posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = params,
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , 1:3], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = params)

color_scheme_set("pink")
mcmc_pairs(posterior, pars = params,
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = params, 
         lags = 50)

# ===========================================================================================
# Age-specific incidence with Vaccine:
# ===========================================================================================
age_group <- c("15-16", "17", "18", "19", "20", "21", "22-23", "24-29", "30-49", "50≥")
incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, .1, .05, .01, .005))

# Convert rate to probability:
p_Age <- 1 - exp(-incidence * 1)

Pr_Age <- cbind(age_group, round(p_Age, 4))
Pr_Age

p_Age <- c(rep(p_Age[1], length(15:16)), rep(p_Age[2], 1), rep(p_Age[3], 1), 
           rep(p_Age[4], 1), rep(p_Age[5], 1), rep(p_Age[6], 1), 
           rep(p_Age[7], length(22:23)), rep(p_Age[8], length(24:29)), 
           rep(p_Age[9], length(30:49)), rep(p_Age[10], length(50:100)))

OddsEff <- mean(OR)
p_AgeVaccinated <- rep(NA, length(p_Age))

for (i in 1:length(p_Age)) {
 p_AgeVaccinated[i] <- ((OddsEff * p_Age[i] / (1 - p_Age[i])) / 
                      (1 + OddsEff * p_Age[i] / (1 - p_Age[i])))
}
# Age dependet probability of infection given vaccine:
p_AgeVaccinated
# Rate of infection:
(-(1 / 1) * log(1 - p_AgeVaccinated))

 # r = - (1 / t) * log(1 - p)



