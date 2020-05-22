library(tidyverse)
library(dampack)
library(readxl)

# ==========================================================================================
# Misc --------------------------------------------------------------------
# ==========================================================================================
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

# ==========================================================================================
# Age-specific all cause mortality --------------------------------------------------
# ==========================================================================================
# Import ASSA mortality table:
mort_data <- read_excel("mortality tables.xls", 
                        sheet = "ASSA data", range = "B3:C94")

# Total Pop. divided by total deaths:
v_r_mort_by_age <- mort_data[, 2] / mort_data[, 1]
v_r_mort_by_age[1:90, ]

# Therefore, all-less HPV mortality calculated as All-cause mortality - HPV mortality, i.e.
# the number of people who die due to cervical cancer.

# Vector of probabilities directly computed with model.

# ===========================================================================================
# Vaccine efficacy --------------------------------------------------------
# ===========================================================================================
# Vaccine efficacy is expressed as the reduction in the occurrence of HPV due to vaccination.
# Individuals who are vaccinated fully experience a rate of occurrence of HPV that is 
# (1 – alpha) times that of those who are not vaccinated (Favato G., et al. 2012).

# COMPARISON GROUPS: PLACEBO (A) AND VACCINE(B)

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
plot(density(rbeta(n = 10000, shape1 = 50.4, shape2 = 12.6)))

# ===========================================================================================
# Vaccine compliance --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:

# Assume full compliance; chi = 1

# ===========================================================================================
# Cross-protection effect --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.
lnorm_params(m = .0740, v = .05)

# ===========================================================================================
# Annual Screening Coverage  --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. HPV and Related Diseases Report
# Proportion of pop. screening coverage every three years according to age group, converted 
# to an annual rate, and then annual probability:

# 12.9% (18-29 years)
- (1 / 3) * log(1 - .129)
1 - exp(-.04603777 * 1)
beta_params(mean = .05, sigma = .05)

# 21.4% (30-39 years)
- (1 / 3) * log(1 - .214)
1 - exp(-.08026616 * 1)
beta_params(mean = .08, sigma = .05)

# 11.5% (40-49 years)
- (1 / 3) * log(1 - .115)
1 - exp(-.04072254 * 1)
beta_params(mean = .04, sigma = .05)

# 7.7% (50-59 years)
- (1 / 3) * log(1 - .077)
1 - exp(-.02670868 * 1)
beta_params(mean = .03,sigma = .05)

# 5.8% (60-69 years)
- (1 / 3) * log(1 - .058)
1 - exp(-.01991667 * 1)
beta_params(mean = .02, sigma = .05)

# ==========================================================================================
# Age-specific Exposure Probability ------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# Anne Bakilana (2005) Age at sexual debut in South Africa, African Journal ofAIDS Research, 
# 4:1, 1-5, DOI: 10.2989/16085900509490335

# Assumes cons. rate.

# From Healthy to Exposed:
p_Sexual_Activity <- rateConversionCons(r = .09, t = 75)
p_Sexual_Activity

plot(p_Sexual_Activity, type = "l", lwd = 2, ylab = "Cum. Pr(Sexual Activity)", 
     xlab = "Ages 15-85", xgap.axis = 100)

# Vector of probabilities directly computed with model.

# ==========================================================================================
# Age-specific Infection Prevalence ----------------------------------------------------
# ==========================================================================================
# Note: will use a relative risk parameter for proportion of HIV+ population who have 
# increased risk of being reinfected rather than clearing HPV. See Italian study.

# Informative prior -------------------------------------
# Mix of the following two sources:
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:

# 2. Data collected from HPV and Related Diseases Report. Used age groups ≥=25.

# Age groups:
age_group <- c("15-16", "17", "18", "19", "20", "≤25", "25-34", 
               "35-44", "45-54", "55≥")
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
# Age-specific Regression of Infection to Exposure -------------------------------
# ==========================================================================================
# Informative prior:

# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# Ages 15-24 years:
beta_params(mean = .719, sigma = .05)
# Ages 25-29 years:
beta_params(mean = .699, sigma = .05)
# Ages 30-39 years:
beta_params(mean = .35, sigma = .05)
# Ages 40-49 years:
beta_params(mean = .201, sigma = .05)
# Ages ≥50 years:
beta_params(mean = .099, sigma = .05)

# ==========================================================================================
# LSIL -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# From Infection progression to Low-Grade Squamous Intraepithelial Lesions:
beta_params(mean = .11, sigma = .1)

# ==========================================================================================
# HSIL -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# From Infection progression to High-Grade Squamous Intraepithelial Lesions:
beta_params(mean = .022, sigma = .01)

# ===========================================================================================
# Subject to No-Treatment  --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1.Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# From LSIL to Clearance:
beta_params(mean = .71, sigma = .05)
# From LSIL to HSIL:
beta_params(mean = .224, sigma = .05)
# From HSIL to Clearance:
beta_params(mean = .355, sigma = .05)
# From HSIL to LSIL:
beta_params(mean = .25, sigma = .05)
# From HSIL to Cancer:
beta_params(mean = .05, sigma = .05)

# ===========================================================================================
# Subject to Treatment (Conization)  --------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# Favato G., et al. 2012: Novel Health Economic Evaluation of a Vaccination Strategy to 
# Prevent HPV-related Diseases.

# From LSIL to Clearance:
beta_params(mean = .8990, sigma = .05)
# From LSIL to HSIL:
beta_params(mean = .12, sigma = .05)
# From HSIL to Clearance:
beta_params(mean = .861, sigma = .05)
# From HSIL to cancer:
beta_params(mean = .015, sigma = .05)

# Probability of Conization
# Immediate:
beta_params(mean = .302, sigma = .05)
# Delayed:
beta_params(mean = .17, sigma = .05)

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

# HIV+ women positive to HPV: 205/277
# HIV- women positive to HPV: 76/207

rHIV_pos <- 205
nHIV_pos <- 277
rHIV_neg <- 76
nHIV_neg <- 207

# Proportion of pop. at increased risk due to HIV:
zeta_alpha <- beta_params(mean = .131, sigma = .05)$alpha
zeta_beta <- beta_params(mean = .131, sigma = .05)$beta

# End file ----------------------------------------------------------------