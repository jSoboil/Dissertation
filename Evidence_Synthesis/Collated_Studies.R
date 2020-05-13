library(tidyverse)
library(readxl)

# ==========================================================================================
# Misc --------------------------------------------------------------------
# ==========================================================================================
# COMPARISON GROUPS: VACCINE(A) AND PLACEBO (B).

# Function to convert rate to probability assuming cons. rate:
rateConversionCons <- function(r, timelength) {
 for (i in 1:timelength) {
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
# Age-specific HPV incidence ----------------------------------------------------
# ==========================================================================================

# Study 1: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Age-specific incidence:
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
barplot(p_Age, names.arg = "From Age 15 to 100", ylab = "Probability of HPV+")

# ==========================================================================================
# Vaccine efficacy --------------------------------------------------------
# ==========================================================================================
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

# Study 1: Villa L, et al. 2006 -------------------------------------------
# High sustained efficacy of a prophylactic quadrivalent human papillomavirus types 
# 6/11/16/18 L1 virus-like particle vaccine through 5 years of follow-up.

# The study enrolled nonpregnant, healthy women who had no prior abnormal Pap smears, and 
# reported a lifetime history of four or fewer male sex partners.

# Events by group:
# tB:0
# nB:235

# tA:3
# nA:233

# Study 2: Garland S., et al. FUTURE I 2007 -------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# Events by group:
# tB:0
# nB:2261

# tA:48
# nA:2279

# Study 3: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tB:0
# nB:2429

# tA:45
# nA:2396

# ==========================================================================================
# LSIL -----------------------------------------------------
# ==========================================================================================
# Low-Grade Squamous Intraepithelial Lesions

# Data collated from per-protocol population. CIN 1 graded as LSIL.

# Study 1: Wei L., et al. 2018 --------------------------------------------
# Efficacy of quadrivalent human papillomavirus vaccine against persistent infection and 
# genital disease in Chinese women: A randomized, placebo-controlled trial with 78-month 
# follow-up.

# Pregnant women and those with a history of genital warts or significant cervical disease, 
# active cervical disease, or prior HPV vaccine recipients were excluded.

# Events by group:
# tB:1
# nB:1271

# tA:29
# nA:1243

# Study 2: Catellsague X., et al. 2011 ------------------------------------
# End-of-study safety, immunogenicity, and efficacy of quadrivalent HPV (types 6, 11, 16, 
# 18) recombinant vaccine in adult women 24–45 years of age.

# Study enrolled 3819 24–45-year-old women with no history of cervical disease or genital 
# warts in the past 5 years.

# Events by group:
# tB:0
# nB:1578

# tA:25
# nA:1583

# Study 3: Garland S., et al. FUTURE I 2007 -------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# Events by group:
# tB:0
# nB:2241

# tA:49
# nA:2258

# Study 4: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tB:2
# nB:2415

# tA:29
# nA:2377

# ==========================================================================================
# HSIL -----------------------------------------------------
# ==========================================================================================
# High-Grade Squamous Intraepithelial Lesions.

# Data collated from per-protocol population. CIN 2/3 graded as HSIL.

# Study 1: Wei L., et al. 2018 --------------------------------------------
# Efficacy of quadrivalent human papillomavirus vaccine against persistent infection and 
# genital disease in Chinese women: A randomized, placebo-controlled trial with 78-month 
# follow-up.

# Pregnant women and those with a history of genital warts or significant cervical disease, 
# active cervical disease, or prior HPV vaccine recipients were excluded.

# Events by group:
# tB:0 
# nB:1271

# tA:3
# nA:1243

# Study 2: Catellsague X., et al. 2011 ------------------------------------
# End-of-study safety, immunogenicity, and efficacy of quadrivalent HPV (types 6, 11, 16, 
# 18) recombinant vaccine in adult women 24–45 years of age.

# Study enrolled 3819 24–45-year-old women with no history of cervical disease or genital 
# warts in the past 5 years.

# Events by group:
# tB:0
# nB:1578

# tA:0
# nA:1583

# Study 3: Garland S., et al. FUTURE I 2007 -------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# Events by group:
# tB:0
# nB:2241

# tA:38
# nA:2258

# Study 4: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tB:1
# nB:2415

# tA:21
# nA:2377

# ==========================================================================================
# Cervical Adenocarcinoma -------------------------------------------------
# ==========================================================================================
# Use observational studies to obtain incidence or prevalence.

# Study 1: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Progression rate HSIL to stage I cancer:

# Data collated from per-protocol population.

# Study 2: Garland S., et al. FUTURE I 2007 -------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# Events by group:
# tB:0
# nB:2241

# tA:6
# nA:2258

# ==========================================================================================
# Median time to response -------------------------------------------------
# ==========================================================================================

