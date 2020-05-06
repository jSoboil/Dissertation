library(tidyverse)

# Function to convert rate to probability assuming cons rate:
rateConversionCons <- function(r, timelength) {
 for (i in 1:timelength) {
  P <- 1 - exp(-r * 1:i)
 }
 print(P)
}

# Each cohort collected from each study's end-point. All data extracted only includes 
# estimates from the per-protocol pop.

# COMPARISON GROUPS: VACCINE(A) AND PLACEBO (B).

# Note on conversion between rates and probabilities (Fleurence R., 2007):
# If an event occurs at a constant rate r per time unit t, then the probability that an 
# event will occur during time t is given by equation 1 (note that the unit of time used
# in r and t must be the same):

#      p = 1 - exp(-r * t)

# On the other hand, if we have a probability and we want to convert it to a rate,

#      r = - (1 / t) * log(1 - p)

# ==========================================================================================
# Age-specific mortality --------------------------------------------------
# ==========================================================================================
# Import WHO lifetable:
mort_data <- read_csv("mort_data.csv", col_types = cols(Female = col_double()), 
    skip = 1)
v_r_mort_by_age <- mort_data[, 2:3]
View(v_r_mort_by_age)

# Age-specific probability of dying when Healthy (all-cause mortality)
p_HDage  <- 1 - exp(-v_r_mort_by_age[1:19, 2] * 1)
p_HDage

# ==========================================================================================
# Probability of Age-specific incidence ----------------------------------------------------
# ==========================================================================================
# Data collated from per-protocol population.

# Study 1: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Age-specific incidence:
age <- c("15-16", "17", "18", "19", "20", "21", "22-23", "24-29", "30-49", "50≥")
incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, .1, .05, .01, .005))

# Convert rate to probability:
p_Age <- 1 - exp(-incidence * 1)
cbind(age, round(p_Age, 4))
barplot(p_Age, names.arg = age)

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
# tA:0
# nA:235

# tB:3
# nB:233

# Study 2: Garland S., et al. 2007 ----------------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# Events by group:
# tA:0
# nA:2261

# tB:48
# nB:2279

# ==========================================================================================
# Probability of LSIL -----------------------------------------------------
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
# tA:1
# nA:1271

# tB:29
# nB:1243

# Study 2: Catellsague X., et al. 2011 ------------------------------------
# End-of-study safety, immunogenicity, and efficacy of quadrivalent HPV (types 6, 11, 16, 
# 18) recombinant vaccine in adult women 24–45 years of age.

# Study enrolled 3819 24–45-year-old women with no history of cervical disease or genital 
# warts in the past 5 years.

# Events by group:
# tA:0
# nA:1578

# tB:25
# nB:1583

# Study 3: Garland S., et al. 2007 ----------------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# Events by group:
# tA:0
# nA:2241

# tB:49
# nB:2258

# ==========================================================================================
# Probability of HSIL -----------------------------------------------------
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
# tA:0 
# nA:1271

# tB:3
# nB:1243

# Study 2: Catellsague X., et al. 2011 ------------------------------------
# End-of-study safety, immunogenicity, and efficacy of quadrivalent HPV (types 6, 11, 16, 
# 18) recombinant vaccine in adult women 24–45 years of age.

# Study enrolled 3819 24–45-year-old women with no history of cervical disease or genital 
# warts in the past 5 years.

# Events by group:
# tA:0
# nA:1578

# tB:0
# nB:1583

# Study 3: Garland S., et al. 2007 ----------------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# Events by group:
# tA:0
# nA:2241

# tB:38
# nB:2258

# ==========================================================================================
# Cervical Adenocarcinoma -------------------------------------------------
# ==========================================================================================
# Data collated from per-protocol population.

# Study 1: Garland S., et al. 2007 ----------------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent Anogenital Diseases.

# Vaccine Efficacy against External Anogenital, Vaginal, and Cervical Lesions Associated 
# with HPV-6, HPV-11, HPV-16, or HPV-18 or Regardless of HPV Type.

# Events by group:
# tA:0
# nA:2241

# tB:6
# nB:2258
















