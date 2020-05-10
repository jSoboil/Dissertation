library(tidyverse)
library(readxl)

# ==========================================================================================
# Thoughts on Model -------------------------------------------------------
# ==========================================================================================



# ==========================================================================================
# Misc --------------------------------------------------------------------
# ==========================================================================================

# Function to convert rate to probability assuming cons rate:
rateConversionCons <- function(r, timelength) {
 for (i in 1:timelength) {
  P <- 1 - exp(-r * 1:i)
 }
 print(P)
}

# Each cohort collected from each study's end-point. All data extracted only includes 
# estimates from per-protocol pop.

# COMPARISON GROUPS: VACCINE(A) AND PLACEBO (B).

# Note on conversion between rates and probabilities (Fleurence R., 2007):
# If an event occurs at a constant rate r per time unit t, then the probability that an 
# event will occur during time t is given by equation 1 (note that the unit of time used
# in r and t must be the same):

#      p = 1 - exp(-r * t)

# On the other hand, if we have a probability and we want to convert it to a rate,

#      r = - (1 / t) * log(1 - p)

# Included statistics for HPV-6, HPV-11, HPV-16, and HPV-18 as model assess quadrivalent
# vaccine efficacy.

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

# Data collated from per-protocol population.

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
length(p_Age)
# I THINK IT IS BEST TO RUN THE VIRTUAL COHORT FROM AGE OF VACCINATION TO 100.
barplot(p_Age, names.arg = "From Age 15 to 100", ylab = "Probability of HPV+")

# Study 2: Richter K., et al. 2013 ----------------------------------------
# Age-specific prevalence of cervical human papillomavirus infection and cytological 
# abnormalities in women in Gauteng Province, South Africa

# This cross-sectional study describes the age-specific prevalence of human papillomavirus 
# (HPV) infection and cytological abnormalities among this urban and peri-urban population.
age_group_2 <- c("<25", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", 
                 "≥55", "Total")
# Total subjects:
N_2 <- c(100, 186, 233, 215, 206, 180, 146, 179, 1445)
# Total HPV positive:
HPV_pos <- c(85, 149, 184, 157, 159, 122, 93, 129, 1084)

cbind(age_group_2, HPV_pos, N_2)

barplot(HPV_pos[-9]/N_2[-9], names.arg = age_group_2[-9], ylab = "Proportion HPV+",
        xlab = "By Age group", main = "Richter K., et al. 2013.", ylim = 0:1)

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

# Study 3: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tA:0
# nA:2429

# tB:45
# nB:2396

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

# Study 4: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tA:2
# nA:2415

# tB:29
# nB:2377

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

# Study 4: Perez G., et al. 2008 ------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus (types 6, 11, 16, 
# 18) L1 virus-like-particle vaccine in Latin American women.

# This study did not exclude subjects with prior HPV infection.

# Events by group:
# tA:1
# nA:2415

# tB:21
# nB:2377

# ==========================================================================================
# Cervical Adenocarcinoma -------------------------------------------------
# ==========================================================================================
# Use observational studies to obtain incidence or prevalence.

# Study 1: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Progression rate HSIL to stage I cancer:

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