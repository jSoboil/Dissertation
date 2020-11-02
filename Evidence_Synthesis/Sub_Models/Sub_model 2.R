library(rjags)
library(R2jags)
library(bayesplot)
library(tidyverse)
library(dampack)
library(MCMCpack)
library(Compositional)

# ==========================================================================================
# Cancer State Progression ------------------------------------------------
# ==========================================================================================

# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.


# From Stage I Cancer to Stage II Cancer
# Probability of progression every 4 years = .9
# Convert to yearly rate
- (1 / 4) * log(1 - .9)
# Yearly probability from Stage I to Stage II:
1 - exp(-0.5756463 * 1)

# Annual probability of symptoms from Stage I to Treatment:
.15

# Yearly probability from Stage I to Stage I:
1 - ((1 - exp(-0.5756463 * 1)) + .15)

# To check probability laws, uncomment the line below:
# 0.4123413 + ((1 - exp(-0.5756463 * 1)) + .15)

# Distribution of Stage I progression:
# alpha = (0.4376587 + 0.15 + 0.4123413)
# Stage_I ~ ddirch(alpha)
alpha.Stage_I <- c(0.4376587, 0.15, 0.4123413)


# From Stage II Cancer to Treatment or Stage III Cancer
# Probability of progression every 3 years = .9
# Convert to yearly rate
- (1 / 3) * log(1 - .9)
# Yearly probability from Stage I to Stage II:
1 - exp(-0.7675284 * 1)

# Annual probability of symptoms from Stage II to Treatment:
0.225

# Yearly probability from Stage II to Stage II:
1 - ((1 - exp(-0.7675284 * 1)) + 0.225)

# To check probability laws, uncomment the line below:
# 0.2391589 + ((1 - exp(-0.7675284 * 1)) + 0.225)

# Distribution of Stage II progression:
# alpha = (0.5358411 + 0.225 + 0.2391589)
# Stage_II ~ ddirch(alpha)
alpha.Stage_II <- c(0.5358411, 0.225, 0.2391589)


# From Stage III Cancer:
# Probability of progression every 2 years = .9
# Convert to yearly rate
- (1 / 2) * log(1 - .9)
# Yearly probability from Stage III to Stage IV:
1 - exp(-1.151293 * 1)

# Yearly probability of symptoms from Stage II to Treatment:
# This doesn't make sense unless it is a proportion of the leftover cohort less those who
# progress to Stage IV...
((1 - (1 - exp(-1.151293 * 1))) * .6)

# Yearly probability from Stage III to Stage III:
1 - (((1 - (1 - exp(-1.151293 * 1))) * .6) + (1 - exp(-1.151293 * 1)))

# To check probability laws, uncomment the line below:
#  0.126491 + (((1 - (1 - exp(-1.151293 * 1))) * .6) + (1 - exp(-1.151293 * 1)))

# Distribution of Stage III progression:
# alpha = (0.6837724 + 0.1897366 + 0.126491)
# Stage_III ~ ddirch(alpha)
alpha.Stage_III <- c(0.6837724, 0.1897366, 0.126491)

# From Stage IV Cancer to Treatment
# Yearly probability of symptoms from Stage IV to Treatment:
.9
 
# Yearly probability from Stage IV to Stage IV:
1 - 0.9
# To check probability laws, uncomment the line below:
# 0.1 + 0.9

# Distribution of Stage IV progression:
# alpha = (0.9 + 0.1)
# Stage_IV ~ ddirch(alpha)
alpha.Stage_IV <- beta_params(mean = .9, sigma = .05)$alpha
beta.Stage_IV <- beta_params(mean = .9, sigma = .05)$beta

# ==========================================================================================
# Vaccine efficacy data ----------------------------------------------------
# ==========================================================================================
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

mu.a.log <- log(Prevalence)

# Optional 'emprical bayes' method:
# mu.a.log <- lnorm_params(m = Prevalence, v = .01)$mu
# sigma.a.log <- lnorm_params(m = Prevalence, v = .01)$sigma
# prec.age <- 1 / (sigma.a.log * sigma.a.log)
# prec.age

# ==========================================================================================
# Sub-model 1 ------------------------------------------------------
# ==========================================================================================
model_String <- "
model {

# SUB-MODEL 1: AGE-SPECIFIC PREVALENCE
# Model parameters abbreviated by .age. Note: this is equivalent to Monte Carlo PSA, 
# as it is technically sampling directly from a prior and is not propogated into a 
# posterior using a likelihood model. However, a hyperprior is used for the population 
# variance to account for a greater uncertainty in each age-population.
  for (i in 1:10) {
    # Monte Carlo:
    omega.age[i] ~ dlnorm(mu.a.log[i], prec.age[i])
    
    # Note in use of pow() function, using -2 is a shorthand inverse
    # method equivalent to 1 / x^2.
    log(prec.age[i]) <- pow(sigma.age[i], -2)
    # Prior on variance for each age group. Note use of half Student-t to draw
    # variance away from 0. See Gelman (2006):
    sigma.age[i] ~ dt(0, eta.age, 1)T(0, )
  }
  
   # Wide hyper-prior on prior variance parameter for SUB-MODEL 1:
   eta.age ~ dunif(0, 1000)
 
# END OF SUB-MODEL 1.

# SUB-MODEL 2: VACCINE-EFFICACY.
# Model parameters abbreviated by .vac.
  for (i in 1:Nstud.vac) {
    # Likelihood:
    rA.vac[i] ~ dbin(pA.vac[i], nA.vac[i])
    rB.vac[i] ~ dbin(pB.vac[i], nB.vac[i])
    
    # Random Effect Logistic model:
    logit(pA.vac[i]) <- mu.vac[i]
    logit(pB.vac[i]) <- mu.vac[i] + delta.vac[i]
    
    # Average effect prior for sub-model 2:
    mu.vac[i] ~ dnorm(0, 1e-6)
    # Prior for sub-model 2 (Random. pop. effect):
    delta.vac[i] ~ dnorm(psi.vac, prec.vac)
  }
  
   # Hyperpriors for SUB-MODEL 2:
   psi.vac ~ dnorm(0, 1.0e-6)
   prec.vac <- pow(tau.vac, -2)
   tau.vac ~  dunif(0, 10)
  
  # Transformations for SUB-MODEL 2:
   # Convert LOR to OR
   OR.vac <- exp(psi.vac)
   # Convert OR to probability
   # for vaccine efficacy
   pEfficacy.vac <- 1 / (1 + OR.vac)

# END OF SUB-MODEL 2.

# SUB-MODEL 3: CANCER PROGRESSION STAGES I-IV.
# Model parameters abbreviated by .canc. Note: this is equivalent to a standard 
# Monte Carlo PSA, as it is technically sampling directly from a prior and it is
# *not* propogated into a posterior using a likelihood model. 

  # Monte Carlo:
   Stage.I.canc ~ ddirch(alpha.Stage_I)
   Stage.II.canc ~ ddirch(alpha.Stage_II)
   Stage.III.canc ~ ddirch(alpha.Stage_III)
   Stage.IV.canc ~ dbeta(alpha.Stage_IV, beta.Stage_IV)

# END OF SUB-MODEL 3.

 }
"
writeLines(text = model_String, con = "CancerStages.txt")

# Transform data into list format so that can be read by JAGS:
data_JAGS <- list(
 # Vaccine efficacy data:
  Nstud.vac = Nstud.vac, 
  rA.vac = rA.vac, nA.vac = nA.vac, 
  rB.vac = rB.vac, nB.vac = nB.vac,
  # Population prevalence:
  mu.a.log = mu.a.log,
  
  # Cancer Stage alpha's:
  alpha.Stage_I = alpha.Stage_I, alpha.Stage_II = alpha.Stage_II,
  alpha.Stage_III = alpha.Stage_III, 
  alpha.Stage_IV = alpha.Stage_IV, beta.Stage_IV = beta.Stage_IV
)

# Parameters to monitor:
params <- c(
  # Vaccine efficacy parameters:
  "OR.vac", "pEfficacy.vac", 
  # Age prevalence:
  "omega.age",
  # Cancer Stage Progression:
  "Stage.I.canc", "Stage.II.canc",
  "Stage.III.canc", "Stage.IV.canc"
  )

# Set no. of iterations, burn-in period and thinned samples:
n.iter <- 25000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 250)

# Run MCMC model:
mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, 
                 model.file = "CancerStages.txt", n.chains = 4, 
                 n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
mod_JAGS

attach.jags(mod_JAGS)

# Visual Inspection of Posterior ------------------------------------------
# Create array of posterior samples:
posterior <- as.array(mod_JAGS$BUGSoutput$sims.array)
dimnames(posterior)

# Inspect trace of chain for each parameter:
color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("pEfficacy.vac", "Stage.I.canc[1]", "Stage.I.canc[2]", 
                               "Stage.I.canc[3]", "Stage.II.canc[1]", "Stage.II.canc[2]", 
                               "Stage.II.canc[3]", "Stage.III.canc[1]", "Stage.III.canc[2]",
                               "Stage.III.canc[3]", "Stage.IV.canc"),
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("pEfficacy.vac", "Stage.I.canc[1]", "Stage.I.canc[2]", 
                                  "Stage.I.canc[3]", "Stage.II.canc[1]", "Stage.II.canc[2]", 
                               "Stage.II.canc[3]", "Stage.III.canc[1]", "Stage.III.canc[2]",
                               "Stage.III.canc[3]", "Stage.IV.canc")
                  )

color_scheme_set("mix-teal-pink")
mcmc_dens(posterior, pars = c("pEfficacy.vac", "Stage.I.canc[1]", "Stage.I.canc[2]", 
                                  "Stage.I.canc[3]", "Stage.II.canc[1]", "Stage.II.canc[2]", 
                               "Stage.II.canc[3]", "Stage.III.canc[1]", "Stage.III.canc[2]",
                               "Stage.III.canc[3]", "Stage.IV.canc[1]", "Stage.IV.canc[2]")
                  )

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("pEfficacy.vac", "Stage.I.canc[1]", "Stage.I.canc[2]", 
                                  "Stage.I.canc[3]", "Stage.II.canc[1]", "Stage.II.canc[2]", 
                               "Stage.II.canc[3]", "Stage.III.canc[1]", "Stage.III.canc[2]",
                               "Stage.III.canc[3]", "Stage.IV.canc[1]", "Stage.IV.canc[2]"),
         lags = 50)

mcmc_combo(posterior, pars = c("pEfficacy.vac", "Stage.III.canc[1]", "Stage.III.canc[2]",
                               "Stage.III.canc[3]", "Stage.IV.canc[1]", "Stage.IV.canc[2]"))

# Add plot for meta analysis results: found in BCEA pg 52:
# ggplot(tr.eff) + 
#  geom_point(aes(y=mean,x=interventions)) + 
#  geom_errorbar(aes(x=interventions,ymin=low,ymax=high),width =0.15) + 
#  coord_flip() + theme_bw() + 
#  labs(x="",y="Probability of smoking cessation",title="Meta-analysis results")

# Cancer State distribution visualisation ------------------------------
diri.contour(a = alpha.Stage_I, x = Stage.I.canc)
diri.contour(a = alpha.Stage_II, x = Stage.II.canc)
diri.contour(a = alpha.Stage_III, x = Stage.III.canc)
plot(density(Stage.IV.canc), lty = 4, lwd = 2, col = "navyblue", xlim = c(.66, 1))

# Example: say 21% of some age pop. group at has HPV, P(No HPV | Vaccination), assuming
# 50% of pop is covered:
(.5 *.21) * apply(pEfficacy.vac, 2, mean)
# ≈ 19.5% protected by vaccine.

# End file ----------------------------------------------------------------