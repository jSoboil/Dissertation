library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)

# The following sub-model simulates age-specific prevalence and average vaccine efficiacy 
# using data from multiple RCTs. Assumed any end-point - infection or disease - signalling 
# efficacy.

# ===========================================================================================
# Collated studies on vaccine efficacy -------------------------------------------------
# ===========================================================================================
# Study 1: Apter et al. (2015)   ------------------------------------------
# Efficacy of Human Papillomavirus 16 and 18 (HPV-16/18) AS04- Adjuvanted Vaccine against 
# Cervical Infection and Pre-cancer in Young Women: Final Event-Driven Analysis of the 
# Randomized, Double- Blind PATRICIA Trial.

# Data collected from 6-month HPV 16/18 ATP cohort.

# Control (A):
# n = 5375
# cases = 435

# Vaccine (B):
# n = 5406
# cases = 32

# Study 2: Harper et al. (2004) -------------------------------------------
# Efficacy of a bivalent L1 virus-like particle vaccine in prevention of infection with human 
# papillomavirus types 16 and 18 in young women: a randomised controlled trial.

# Data collected from 27 month HPV 16/18 ATP cohort.

# Control (A):
# n = 355
# cases = 41

# Vaccine (B):
# n = 366
# cases = 12

# Study 3: Harper et al. (2006) -------------------------------------------
# Sustained efficacy up to 4·5 years of a bivalent L1 virus-like particle vaccine against 
# human papillomavirus types 16 and 18: follow-up from a randomised control trial.

# Data collected from 54 month HPV 16/18 ATP cohort.

# Control (A):
# n = 277
# cases = 28

# Vaccine (B):
# n = 310
# cases = 1

# Study 4: Herrero et al. (2011) ------------------------------------------
# Prevention of Persistent Human Papillomavirus Infection by an HPV16/18 Vaccine: A 
# Community-Based Randomized Clinical Trial in Guanacaste, Costa Rica.

# Data collected from 22-34 month follow-up HPV 16/18 ATP cohort.

# Control (A):
# n = 2239
# cases = 38

# Vaccine (B):
# n = 2190
# cases = 3

# Study 5: Herrero et al. (2013) ------------------------------------------
# Reduced Prevalence of Oral Human Papillomavirus (HPV) 4 Years after Bivalent HPV Vaccination
# in a Randomized Clinical Trial in Costa Rica.

# Data collected from 48 month 16/18 group. Note!: likely need to discount this study study has
# combined ATP and other cohorts.

# Control (A):
# n = 2924
# cases = 219

# Vaccine (B):
# n = 2910
# cases = 61

# Study 6: Konno et al. (2010) --------------------------------------------
# Efficacy of Human Papillomavirus Type 16/18 AS04YAdjuvanted Vaccine in Japanese Women Aged 
# 20 to 25 Years.

# Data collected from 6 month HPV 16/18 ATP cohort.

# Control (A):
# n = 392
# cases = 15

# Vaccine (B):
# n = 387
# cases = 0

# Study 7: Naud et al. (2014) ---------------------------------------------
# Sustained efficacy, immunogenicity, and safety of the HPV-16/18 AS04-adjuvanted vaccine: 
# Final analysis of a long-term follow-up study up to 9.4 years post-vaccination.

# Data collected from 12 month persistent infection HPV 16/18 ATP cohort, using the combined
# analysis from all follow-up studies, HPV-001/007/023.

# Control (A):
# n = 175
# cases = 10

# Vaccine (B):
# n = 193
# cases = 0

# Study 8: Paavonen et al. (2007) -----------------------------------------
# Efficacy of a prophylactic adjuvanted bivalent L1 virus-like-particle vaccine against 
# infection with human papillomavirus types 16 and 18 in young women: an interim analysis of a
# phase III double-blind, randomised controlled trial.

# Data collected from primary endpoint CIN2+ HPV 16/18.

# Control (A):
# n = 7838
# cases = 21

# Vaccine (B):
# n = 7788
# cases = 2

# Study 9: Paavonen et al. (2009) -----------------------------------------
# Efficacy of human papillomavirus (HPV)-16/18 AS04-adjuvanted vaccine against cervical 
# infection and precancer caused by oncogenic HPV types (PATRICIA): final analysis of a 
# double-blind, randomised study in young women.

# Data collected from 34·9 month HPB 16/18 ATP cohort. Indicator used CIN2+.

# Control (A):
# n = 7312
# cases = 56

# Vaccine (B):
# n = 7344
# cases = 4

# Study 10: Romanowski et al. (2009) --------------------------------------
# Sustained efficacy and immunogenicity of the human papillomavirus (HPV)-16/18 AS04-adjuvanted
# vaccine: analysis of a randomised placebo-controlled trial up to 6·4 years.

# Data collected from 12-month persistent infection HPV 16/18 ATP cohort.

# Control (A):
# n = 372
# cases = 20

# Vaccine (B):
# n = 401
# cases = 0

# Study 11: Zhu et al. (2014) ---------------------------------------------
# Efficacy, immunogenicity and safety of the HPV-16/18 AS04- adjuvanted vaccine in healthy 
# Chinese women aged 18–25 years: Results from a randomized controlled trial.

# Data collected from 6-Month PI and/or CIN1+ (primary endpoint) HPV 16/18 ATP cohort.

# Control (A):
# n = 2502
# cases = 17

# Vaccine (B):
# n = 2497
# cases = 1

# Vaccine efficacy data ----------------------------------------------------
# PLACEBO (A) AND VACCINE(B).
Nstud.vac <- 11
rA.vac <- c(435, 41, 28, 38, 219, 15, 10, 21, 56, 20, 17)
nA.vac <- c(5375, 355, 277, 2239, 2924, 392, 175, 7838, 7312, 372, 2502)
rB.vac <- c(32, 12, 1, 3, 61, 0, 0, 2, 4, 0, 1)
nB.vac <- c(5406, 366, 310, 2190, 2910, 387, 193, 7788, 7344, 401, 2497)
# Check:
length(rA.vac) == length(nA.vac)
length(rB.vac) == length(nB.vac)



# ==========================================================================================
# Age-specific infection ----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:

# Informative prior sampling model:
# omega.age[i] ~ dlnorm(mu.log[i], prec.log[i])

# Age groups:
age_group <- c("15-16", "17", "18", "19", "20", "≤21", "22-23", 
               "24-29", "30-49", "≥55")

# Estimated Prevalence:
Prevalence <- c(.09516258, .1130796, .139292, .1563352, .139292, 
                .1130796, .09516258, .04877058, .009950166, .004987521)

mu.a.log <-log(Prevalence)

# mu.a.log <- lnorm_params(m = Prevalence, v = sd(Prevalence)^2)$mu
# Optional non-hierarchical model on variance:
# sigma.a.log <- lnorm_params(m = Prevalence, v = sd(Prevalence)^2)$sigma
# prec.age <- 1 / (sigma.a.log * sigma.a.log)
# prec.age

# ==========================================================================================
# Sub-model 1 ------------------------------------------------------
# ==========================================================================================
model_String <- "
model {

# SUB-MODEL 1: AGE-SPECIFIC PREVALENCE,
  # model parameters abbreviated by .age
    # Note: equivalent to Monte Carlo PSA, as it 
    # is technically sampling directly from a prior
    # without any likelihood model.
  for (i in 1:10) {
    omega.age[i] ~ dlnorm(mu.a.log[i], prec.age[i])
    
    # Sigma precision:
    log(prec.age[i]) <- 1 / (sigma.age[i] ^ 2)
    # Sigma Pior:
    sigma.age[i] ~ dt(0, sigma.hyperprior, 1)T(0, )
    
  }
    # Hyper-prior for sigma.age prior:
    sigma.hyperprior ~ dunif(0, 100)

# END OF SUB-MODEL 1.

# SUB-MODEL 2: VACCINE-EFFICACY,
  # model parameters abbreviated by .vac
  for (i in 1:Nstud.vac) {
    # Likelihood:
    rA.vac[i] ~ dbin(pA.vac[i], nA.vac[i])
    rB.vac[i] ~ dbin(pB.vac[i], nB.vac[i])
    # Logistic function:
    logit(pA.vac[i]) <- mu.vac[i]
    logit(pB.vac[i]) <- mu.vac[i] + delta.vac[i]
    
    # Average effect prior for sub-model 2:
    mu.vac[i] ~ dnorm(0, 1e-4)
    # Prior for sub-model 2 (Random. pop. effect):
    delta.vac[i] ~ dnorm(psi.vac, prec.vac)
  }
  
   # Hyperpriors for sub-model 2:
   psi.vac ~ dnorm(0, 1.0e-4)
   prec.vac <- 1 / pow(tau.vac, 2)
   tau.vac ~ dunif(0, 100)
  
  # Transformations for Sub-model 2:
   # Convert LOR to OR
   OR.vac <- exp(psi.vac)
   # Convert OR to probability
   # for vaccine efficacy
   pEfficacy.vac <- 1 / (1 + OR.vac)

# END OF SUB-MODEL 2.

}
"
writeLines(text = model_String, con = "Age_and_Efficacy.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(
 # Vaccine efficacy data:
  Nstud.vac = Nstud.vac, 
  rA.vac = rA.vac, nA.vac = nA.vac, 
  rB.vac = rB.vac, nB.vac = nB.vac,
  # Population prevalence:
  mu.a.log = mu.a.log
)

# Parameters to monitor:
params <- c(
  # Vaccine efficacy parameters:
  "OR.vac", "pEfficacy.vac", 
  # Age prevalence:
  "omega.age"
  )

# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 40000
n.burnin <- 15000
n.thin <- floor((n.iter - n.burnin) / 250)

# Run MCMC model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "Age_and_Efficacy.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS

# Attach JAGS model to global envir:
attach.jags(mod_JAGS)

# Visual Inspection of Posterior ------------------------------------------
# Create array of posterior samples:
posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

# Inspect trace of chain for each parameter:
color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("OR.vac", "pEfficacy.vac", "omega.age[1]", "omega.age[2]"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[, , ], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("pEfficacy.vac", "OR.vac", 
                                      "omega.age[1]", "omega.age[2]", "omega.age[6]")
                  )

color_scheme_set("mix-teal-pink")
mcmc_dens(posterior, pars = c("pEfficacy.vac", "OR.vac", 
                              "omega.age[1]", "omega.age[2]", "omega.age[6]")
                  )

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("pEfficacy.vac", "OR.vac", 
                             "omega.age[1]", "omega.age[2]"),
         lags = 50)
mcmc_acf_bar(posterior, pars = c("pEfficacy.vac", "OR.vac", 
                             "omega.age[1]", "omega.age[2]"),
         lags = 25)

mcmc_areas(posterior, pars = c("pEfficacy.vac","omega.age[2]"))

mcmc_combo(posterior, pars = c("pEfficacy.vac", "OR.vac", 
                             "omega.age[1]", "omega.age[2]"))

# Add plot for meta analysis results: found in BCEA pg 52:
# ggplot(tr.eff) + 
#  geom_point(aes(y=mean,x=interventions)) + 
#  geom_errorbar(aes(x=interventions,ymin=low,ymax=high),width =0.15) + 
#  coord_flip() + theme_bw() + 
#  labs(x="",y="Probability of smoking cessation",title="Meta-analysis results")

# Calibration: ------------------------------------------------------------
# Age-specific HPV prevalence:
plot(apply(omega.age, 2, mean), type = "b", ylab = "Prevalence", xlab = "By Age",
     col = "grey", lwd = 3)

# Example: say 21% of some age pop. group at has HPV, P(No HPV | Vaccination):
.21 * apply(pEfficacy.vac, 2, mean)
# ≈ 19.5% protected by vaccine.

# End file ----------------------------------------------------------------