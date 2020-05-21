library(tidyverse)
library(dampack)
library(readxl)

# ==========================================================================================
# Misc --------------------------------------------------------------------
# ==========================================================================================
# COMPARISON GROUPS: VACCINE(B) AND PLACEBO (A).

# Function to convert rate to probability assuming cons. rate over time:
rateConversionCons <- function(r, t) {
 for (i in 1:t) {
  P <- 1 - exp(-r * 1:i)
 }
 print(P)
}

# Each cohort collected from each study's end-point. All data extracted only includes 
# estimates from per-protocol pop and susceptible population. This is due to the fact that
# the per protocol analysis strategy may be subject to bias as the reasons for non-compliance
# may be related to treatment (Cochrane.org/glossary)

# Note on conversion between rates and probabilities (Fleurence R., 2007):
# If an event occurs at a constant rate r per time unit t, then the probability that an 
# event will occur during time t is given by equation 1 (note that the unit of time used
# in r and t must be the same):

#      p = 1 - exp(-r * t)

# On the other hand, if we have a probability and we want to convert it to a rate,

#      r = - (1 / t) * log(1 - p)

# Note: the included statistics for HPV-6, HPV-11, HPV-16, and HPV-18 as model assesses
# quadrivalent vaccine efficacy.

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
# tA:35
# nA:233
# tB:4
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

# Study 6: Garland S., et al 2007 -----------------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# A phase 3 trial was conducted to evaluate the efficacy of a prophylactic quadrivalent 
# vaccine in preventing anogenital diseases associated with human papillomavirus (HPV) types
# 6, 11, 16, and 18.

# Events by group:
# tA:60
# nA:2279
# tB:0
# nB:2261

# Study 7: Koutsky L., et al. 2007 ----------------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent High-Grade Cervical Lesions

# Randomized, double-blind trial, we assigned 12,167 women between the ages of 15 and 26 
# years to receive three doses of either HPV-6/11/16/18 vaccine or placebo, administered at 
# day 1, month 2, and month 6.

# Events by group:
# tA:42
# nA:5260
# tB:1
# nB:5305

# Study 8: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 
# 16, 18) L1 virus-like-particle vaccine in Latin American women.

# Events by group:
# tA:41
# nA:2377
# tB:3
# nB:2415

# ===========================================================================================
# Vaccine coverage --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:
beta_params(mean = .8, sigma = .05)
# Assumes 80% coverage
plot(density(rbeta(n = 1000, shape1 = 50.4, shape2 = 12.6)))

# ===========================================================================================
# Vaccine compliance --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:
beta_params(mean = .95, sigma = .05)
# Assume full compliance.

# ===========================================================================================
# Efficacy decrease due to non-compliance ------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Favato G., et al. 2012. Bayesian HPV model:
beta_params(mean = .5040, sigma = .05)
# Assume full compliance.

# ===========================================================================================
# Cross-protection effect --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Favato G., et al. 2012. Bayesian HPV model.
lnorm_params(m = .0740, v = .05)

# ==========================================================================================
# Age-specific all cause mortality --------------------------------------------------
# ==========================================================================================
# Import ASSA mortality table:
mort_data <- read_excel("Evidence_Synthesis/mortality tables.xls", 
                        sheet = "ASSA data", range = "B3:C94")

# Total Pop. divided by total deaths:
v_r_mort_by_age <- mort_data[, 2] / mort_data[, 1]
v_r_mort_by_age

# Therefore, all-less HPV mortality calculated as All-cause mortality - HPV mortality, i.e.
# the number of people who die due to cervical cancer.

# ==========================================================================================
# Age-specific infection ----------------------------------------------------
# ==========================================================================================
# Note: will use a relative risk parameter for proportion of HIV+ population who have 
# increased risk of being reinfected rather than clearing HPV. See Italian study.

# Informative prior -------------------------------------
# Mix of the following two sources:
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:

# 2. Data collected from HPV and Related Diseases Report. Used age groups 25+.

# Age groups:
age_group <- c("15-16", "17", "18", "19", "20", "≤25", "25-34", 
               "35-44", "45-54", "55-64")
# Estimated Prevalence:
1 - exp(-.17)
Prevalence <- c(0.09516258, 0.1130796, 0.139292, 0.1563352, 0.139292, 
                .435, .3643, 0.2051, 0.1852, 0.1654)

mu.a.log <- lnorm_params(m = Prevalence, v = .1)$mu
sigma.a.log <- 1 / (lnorm_params(m = Prevalence, v = .1)$sigma) ^ 2
mu.a.log
sigma.a.log
cbind(age_group, Prevalence)
# Informative prior sampling model:
# omega.age[i] ~ dlnorm(mu.log[i], sigma.log[i])

# ==========================================================================================
# Condyloma (Genital Warts) -----------------------------------------------------------
# ==========================================================================================
# Probability of Condyloma.

# Informative prior: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa.

# ==========================================================================================
# LSIL -----------------------------------------------------
# ==========================================================================================
# Low-Grade Squamous Intraepithelial Lesions

# Informative prior: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa.

# Regression probability LSIL to Cleared:
beta_params(mean = .9, sigma = .01)
plot(density(rbeta(n = 100000, shape1 = 809.1, 89.9)))

# ==========================================================================================
# HSIL -----------------------------------------------------
# ==========================================================================================
# High-Grade Squamous Intraepithelial Lesions.

# Informative prior: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa.

# Regression probability from HSIL to LSIL:
# 1 year transition probability (Miller and Homan, 1994)
1 - exp(-.35 / 6)
# Beta parameters:
beta_params(mean = 0.05666455, sigma = .01)
plot(density(rbeta(n = 10000, shape1 = 30.23262, shape2 = 503.3042)))

# Progression probability from HSIL to Cancer:
# 1 year transition probability (Miller and Homan, 1994)
1 - exp(-.4 / 12)
# Beta parameters:
beta_params(mean = 0.0327839, sigma = .01)
plot(density(rbeta(n = 10000, shape1 = 10.3627, shape2 = 305.7285)))

# ==========================================================================================
# Cervical Adenocarcinoma -------------------------------------------------
# ==========================================================================================

# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# Proportion of inviduals in each of the four FIGO stages who progress to cancer, in one
# cycle:
# FIGO I 
# dDirc(.5500)
# FIGO II
# dDirc(.1500)
# FIGO III
# dDirc(.1200)
# FIGO IV
# dDirc(.1800)
alpha.c <- list(.55, .15, .12, .18)

# Survival probabilities over four years, according to each FIGO stage:
# Year 1 survival by stage:
# FIGO I
alphaI.year_1 <- beta_params(mean = .9770, sigma = .05)$alpha
betaI.year_1 <- beta_params(mean = .9770, sigma = .05)$beta
# FIGO II
alphaII.year_1 <- beta_params(mean = .8290, sigma = .05)$alpha
betaII.year_1 <- beta_params(mean = .8290, sigma = .05)$beta
# FIGO III
alphaIII.year_1 <- beta_params(mean = .59, sigma = .05)$alpha
betaIII.year_1 <- beta_params(mean = .59, sigma = .05)$beta
# FIGO IV
alphaIV.year_1 <- beta_params(mean = .5020, sigma = .05)$alpha
betaIV.year_1 <- beta_params(mean = .5020, sigma = .05)$beta

# Year 2 survival by stage:
# FIGO I
alphaI.year_2 <- beta_params(mean = .9790, sigma = .05)$alpha
betaI.year_2 <- beta_params(mean = .9790, sigma = .05)$beta
# FIGO II
alphaII.year_2 <- beta_params(mean = .8330, sigma = .05)$alpha
betaII.year_2 <- beta_params(mean = .8330, sigma = .05)$beta
# FIGO III
alphaIII.year_2 <- beta_params(mean = .6930, sigma = .05)$alpha
betaIII.year_2 <- beta_params(mean = .6930, sigma = .05)$beta
# FIGO IV
alphaIV.year_2 <- beta_params(mean = .7820, sigma = .05)$alpha
betaIV.year_2 <- beta_params(mean = .7820, sigma = .05)$beta

# Year 3 survival by stage:
# FIGO I
alphaI.year_3 <- beta_params(mean = .9630, sigma = .05)$alpha
betaI.year_3 <- beta_params(mean = .9630, sigma = .05)$beta
# FIGO II
alphaII.year_3 <- beta_params(mean = .7550, sigma = .05)$alpha
betaII.year_3 <- beta_params(mean = .7550, sigma = .05)$beta
# FIGO III
alphaIII.year_3 <- beta_params(mean = .7780, sigma = .05)$alpha
betaIII.year_3 <- beta_params(mean = .7780, sigma = .05)$beta
# FIGO IV
alphaIV.year_3 <- beta_params(mean = .7220, sigma = .05)$alpha
betaIV.year_3 <- beta_params(mean = .7220, sigma = .05)$beta

# Year 4 survival by stage:
# FIGO I
alphaI.year_4 <- beta_params(mean = .9890, sigma = .05)$alpha
betaI.year_4 <- beta_params(mean = .9890, sigma = .05)$beta
# FIGO II
alphaII.year_4 <- beta_params(mean = .8690, sigma = .05)$alpha
betaII.year_4 <- beta_params(mean = .8690, sigma = .05)$beta
# FIGO III
alphaIII.year_4 <- beta_params(mean = .9290, sigma = .05)$alpha
betaIII.year_4 <- beta_params(mean = .9290, sigma = .05)$beta
# FIGO IV
alphaIV.year_4 <- beta_params(mean = .9250, sigma = .05)$alpha
betaIV.year_4 <- beta_params(mean = .9250, sigma = .05)$beta

# ==========================================================================================
# Risk Increase (HIV+) ------------------------------------------------------
# ==========================================================================================

# Study: Mbulawa Z., et al. 2015 ------------------------------------------
# Human papillomavirus prevalence in South African women and men according to age and human 
# immunodeficiency virus status.

# Then use proportion of HIV+ pop. for how many at risk at time j. See Italian study.

# Participants were 208 HIV-negative women, 278 HIV-positive women, 325 HIV-negative men and
# 161 HIV-positive men between the ages of 18–66 years. HPV types were determined in cervical
# and penile cells by Roche Linear Array HPV genotyping assay.


# End file ----------------------------------------------------------------