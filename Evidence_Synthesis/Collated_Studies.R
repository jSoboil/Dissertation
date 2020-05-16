library(tidyverse)
library(readxl)

# ==========================================================================================
# Misc --------------------------------------------------------------------
# ==========================================================================================
# COMPARISON GROUPS: VACCINE(A) AND PLACEBO (B).

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

# Informative prior: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Age-specific incidence:
age_group <- c("15-16", "17", "18", "19", "20", "21", "22-23", "24-29", "30-49", "50≥")
incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, .1, .05, .01, .005))
cbind(age_group, incidence)

# Study: Peto J., et al 2004 --------------------------------------------
# Cervical HPV infection and neoplasia in a large population-based prospective study: the 
# Manchester cohort.

# 15-19:
# n = 319
# r = 69

# 20-24:
# n = 466
# r = 92

# 25-29:
# n = 586
# r = 93

# 30-34:
# n = 1191
# r = 86

# 35-39:
# n = 914
# r = 47

# 40-44:
# n = 926
# r = 28

# 45-49:
# n = 508
# r = 15

# 50-54:
# n = 1104
# r = 29

# 55-59:
# n = 114
# r = 1

# ==========================================================================================
# Vaccine efficacy --------------------------------------------------------
# ==========================================================================================
# Vaccine efficacy, expressed as the reduction in the occurrence of HPV due to vaccination.
# Individuals who are vaccinated fully experience a rate of occurrence of HPV that is 
# (1 – alpha) times that of those who are not vaccinated (Favato G., et al. 2012).

# All data collected from quadrivalent vaccine RCTs per protocol pop.

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

# Unrestricted susceptible pop:
# Events by group:
# tB:4
# nB:2667

# tA:81
# nA:2684

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

# Unrestricted susceptible pop:
# Events by group:
# tB:3
# nB:5865

# tA:62
# nA:5863

# Study 7: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle 2vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tB:3
# nB:1990

# tA:25
# nA:1880
	
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