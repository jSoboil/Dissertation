library(tidyverse)
library(dampack)
library(readxl)

# ==========================================================================================
# Misc --------------------------------------------------------------------
# ==========================================================================================
# COMPARISON GROUPS: VACCINE(B) AND PLACEBO (A).

# Function to convert rate to probability assuming cons. rate by time:
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
mort_data <- read_excel("Evidence_Synthesis/mortality tables.xls", 
                        sheet = "ASSA data", range = "B3:C94")

# Total Pop. divided by total deaths:
v_r_mort_by_age <- mort_data[, 2] / mort_data[, 1]
v_r_mort_by_age

# Therefore, all-less HPV mortality calculated as All-cause mortality - HPV mortality, i.e.
# the number of people who die due to cervical cancer.

# ==========================================================================================
# Age-specific HPV infection ----------------------------------------------------
# ==========================================================================================
# Note: will use a relative risk parameter for proportion of HIV+ population who have 
# increased risk of being reinfected rather than clearing HPV. See Italian study.

# Informative prior -------------------------------------
# Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa:
age_group <- c("15-16", "17", "18", "19", "20", "21", "22-23", "24-29", "30-49", "50≥")
incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, .1, .05, .01, .005))

# Data collected from HPV and Related Diseases Report:
age_group <- c("≤25", "25-34", "35-44", "45-54", "55-64")
# Note: probability converted to yearly rate for use in lognormal distribution.
# r = - (1 / t) * log(1 - p)
incidence <- c(.44, .37, .205, .18, .17)
cbind(age_group, incidence)

mu.log <- lnorm_params(m = incidence, v = 1)$mu
sigma.log <- lnorm_params(m = incidence, v = 1)$sigma

# Sampling Model :
# omega.age[i] ~ dlnorm(mu.log[i], sigma.log[i])
plot(density(rlnorm(n = 10000, meanlog = mu.log[4], sdlog = sigma.log[4])))

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

# ==========================================================================================
# Genital Warts -----------------------------------------------------------
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

# Regression probability LSIL to HPV or Cleared:
rateConversionCons(r = .65, t = 6)
1 - exp(-.65 * 6) # 15-34 years
1 - exp(-.4 * 6) # 15-34 years
# Of the regressed, LSIL reverting to Cleared (90%):
.9 * (1 - exp(-.65 * 6))
.9 * (1 - exp(-.4 * 6))
# Progression probability LSIL to HSIL:
1 - exp(-.1 * 6) # 15-34 years.
1 - exp(-.35 * 6) # ≥ 35 years

# ==========================================================================================
# HSIL -----------------------------------------------------
# ==========================================================================================
# High-Grade Squamous Intraepithelial Lesions.

# Informative prior: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa.

# Regression HSIL to LSIL or Cleared:
1 - exp(-.35 * 6)
# Proportion of regressed, HSIL to Cleared:
.5 * 1 - exp(-.35 * 6)
# Progression probability for HSIL to Stage I Cancer:
1 - exp(-.4 * 12)

# ==========================================================================================
# Cervical Adenocarcinoma -------------------------------------------------
# ==========================================================================================

# Informative Prior: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa.

# ==========================================================================================
# Relative Risk HIV+ ------------------------------------------------------
# ==========================================================================================

# Study: Mbulawa Z., et al. 2015 ------------------------------------------
# Human papillomavirus prevalence in South African women and men according to age and human 
# immunodeficiency virus status.

# Then use proportion of HIV+ pop. for how many at risk at time j. See Italian study.

# Participants were 208 HIV-negative women, 278 HIV-positive women, 325 HIV-negative men and
# 161 HIV-positive men between the ages of 18–66 years. HPV types were determined in cervical
# and penile cells by Roche Linear Array HPV genotyping assay.


# End file ----------------------------------------------------------------