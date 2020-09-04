library(tidyverse)
library(dampack)
library(readxl)

# ==========================================================================================
# Misc --------------------------------------------------------------------
# ==========================================================================================
# Function to convert rate to probability assuming cons. rate over time:
rateConv.consR <- function(r, t) {
 for (i in 1:t) {
  P <- 1 - exp(-r * 1:i)
 }
 print(P)
}

# Each cohort collected from each study's end-point. All data extracted only includes 
# estimates from per-protocol pop. However, the per protocol analysis strategy may be 
# subject to bias as the reasons for non-compliance may be related to treatment 
# (Cochrane.org/glossary).

# Note on conversion between rates and probabilities (Fleurence R., 2007):
# If an event occurs at a constant rate r per time unit t, then the probability that an 
# event will occur during time t is given by equation 1 (note that the unit of time used
# in r and t must be the same):

#    p = 1 - exp(-r * t)

# On the other hand, if we have a probability and we want to convert it to a rate,

#    r = - (1 / t) * log(1 - p)

# Note: the included statistics are for HPV-16 and HPV-18 as model assesses bivalent vaccine 
# efficacy.

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

# Study 1: Apter et al. (2015)   ------------------------------------------
# Efficacy of Human Papillomavirus 16 and 18 (HPV-16/18) AS04- Adjuvanted Vaccine against 
# Cervical Infection and Precancer in Young Women: Final Event-Driven Analysis of the 
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

# ===========================================================================================
# Vaccine coverage --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:
beta_params(mean = .8, sigma = .05)
# Assumes 80% coverage with an se of 5%.

# ===========================================================================================
# Vaccine compliance --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
# Used age groups 15-20:

# Assume full compliance; chi = 1

# ==========================================================================================
# Age-specific Infection Prevalence ----------------------------------------------------
# ==========================================================================================
# Note: will use a relative risk parameter for proportion of HIV+ population who have 
# increased risk of being reinfected rather than clearing HPV. See Italian study.

# Informative prior -------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Age groups:
age_group <- c("15-16", "17", "18", "19", "20", "≤21", "22-23", 
               "24-29", "30-49", "≥55")
# Estimated Prevalence:
1 - exp(-.005)
Prevalence <- c(.09516258, .1130796, .139292, .1563352, .139292, 
                .1130796, .09516258, .04877058, .009950166, .004987521)
length(Prevalence)
mu.a.log <- lnorm_params(m = Prevalence, v = .025)$mu
sigma.a.log <- lnorm_params(m = Prevalence, v = .025)$sigma
mu.a.log
sigma.a.log
cbind(age_group, Prevalence)
# Informative prior sampling model:
# omega.age[i] ~ dlnorm(mu.log[i], sigma.log[i])

# ==========================================================================================
# Age-specific Regression of Infection to Normal -------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Ages 15-24 years:
1 - exp(-.7 * 1.5)
beta_params(mean = .65, sigma = .05)
# Ages 25-29 years:
1 - exp(-5 * 1.5)
beta_params(mean = 0.9994, sigma = .05)
# Ages ≥30 years:
1 - exp(-.15 * 1.5)
beta_params(mean = 0.2014838, sigma = .05)

# ==========================================================================================
# From Infection to LSIL -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

1 - exp(-.2 * 3)
beta_params(mean = 0.4511884, sigma = .05)

# ==========================================================================================
# From Infection to HSIL -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

beta_params(mean = .1, sigma = .05)

# ==========================================================================================
# From LSIL to Infection or Normal ----------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Ages 15-34:
1 - exp(-.65 * 6)
beta_params(mean = 0.9797581, sigma = .05)
# Proportion reverting to Normal is .9:
1 - ((0.9797581 - 0.9797581 * .9) + 0.9797581 * .9)

# Ages ≥ 35:
1 - exp(-.4 * 6)
beta_params(mean = 0.909282, sigma = .05)
# Proportion reverting to Normal is .9:
1 - ((0.909282 - 0.909282 * .9) + 0.909282 * .9)

# ==========================================================================================
# From LSIL to HSIL -------------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Ages 15-34:
1 - exp(-.1 * 6)
beta_params(mean = 0.4511884, sigma = .05)

# Ages ≥ 35:
1 - exp(-.35 * 6)
beta_params(mean = 0.8775436, sigma = .05)

# ==========================================================================================
# From HSIL to LSIL or Normal ---------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

1 - exp(-.35 * 6)
# Proportion reverting to Normal is .5:
1 - ((0.8775436 - 0.8775436 * .5) + 0.8775436 * .5)

# ==========================================================================================
# From HSIL to Stage I Cancer -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
1 - exp(-.4 * 10)

# ==========================================================================================
# From Stage I Cancer to Treatment or Stage II Cancer ---------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Probability of symptoms every 4 years = .9
# Convert to yearly rate
- (1 / 4) * log(1 - .9)
# Convert to yearly transition probability to treatment
1 - exp(-.5756463 * 1)

# Yearly probability of progression to Stage II Cancer:
1 - exp(-.15 * 1)

# ==========================================================================================
# From Stage II Cancer to Treatment or Stage III Cancer -------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Probability of symptoms and therefore treatment every 4 years = .9
# Convert to yearly rate
- (1 / 3) * log(1 - .9)
# Convert to yearly transition probability to treatment
1 - exp(-.7675284 * 1)

# Yearly probability of progression to Stage III Cancer:
1 - exp(-.225 * 1)

# ==========================================================================================
# From Stage III Cancer to Treatment or Stage IV Cancer -------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Probability of symptoms and therefore treatment every 4 years = .9
# Convert to yearly rate
- (1 / 2) * log(1 - .9)
# Convert to yearly transition probability to treatment
1 - exp(-1.151293 * 1)

# Yearly probability of progression to Stage IV Cancer:
1 - exp(-.6 * 1)

# ==========================================================================================
# From Stage IV Cancer to Treatment -------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Probability of symptoms and therefore treatment every 4 years = .9
# Convert to yearly rate
- (1 / 2) * log(1 - .9)
# Convert to yearly transition probability to treatment
1 - exp(-1.151293 * 1)

# ==========================================================================================
# From Stage I to Survival for Years 1-5 -----------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Year 1
.9688
beta_params(mean = .9688, sigma = .05)
# Year 2
.9525
beta_params(mean = .9525, sigma = .05)
# Year 3
.9544
beta_params(mean = .9544, sigma = .05)
# Year 4
.9760
beta_params(mean = .9760, sigma = .05)
# Year 5
.9761
beta_params(mean = .9761, sigma = .05)

# ==========================================================================================
# From Stage II to Survival for Years 1-5 -----------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Year 1
.9066
beta_params(mean = .9066, sigma = .05)
# Year 2
.8760
beta_params(mean = .8760, sigma = .05)
# Year 3
.9225
beta_params(mean = .9225, sigma = .05)
# Year 4
.9332
beta_params(mean = .9332, sigma = .05)
# Year 5
.9604
beta_params(mean = .9604, sigma = .05)

# ==========================================================================================
# From Stage III to Survival for Years 1-5 -----------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Year 1
.7064
beta_params(mean = .7064, sigma = .05)
# Year 2
.7378
beta_params(mean = .7378, sigma = .05)
# Year 3
.8610
beta_params(mean = .8610, sigma = .05)
# Year 4
.9231
beta_params(mean = .9231, sigma = .05)
# Year 5
.9142
beta_params(mean = .9142, sigma = .05)

# ==========================================================================================
# From Stage IV to Survival for Years 1-5 -----------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Year 1
.3986
beta_params(mean = .3986, sigma = .05)
# Year 2
.4982
beta_params(mean = .4982, sigma = .05)
# Year 3
.7638
beta_params(mean = .7638, sigma = .05)
# Year 4
.8652
beta_params(mean = .8652, sigma = .05)
# Year 5
.8592
beta_params(mean = .8592, sigma = .05)

# ===========================================================================================
# Annual Cancer Screening Coverage  --------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. HPV and Related Diseases Report
# Proportion of pop. screening coverage every three years according to age group, converted 
# to an annual rate, and then annual probability:

# 18-29 years = 12.9%
- (1 / 3) * log(1 - .129) # 3-year Probability converted to 1-year Rate.
1 - exp(-.04603777 * 1) # 1-year Probability
beta_params(mean = .04499411, sigma = .05)

# 30-39 years = 21.4%
- (1 / 3) * log(1 - .214) # 3-year Probability converted to 1-year Rate.
1 - exp(-.08026616 * 1) # 1-year Probability
beta_params(mean = .07712932, sigma = .05)

# 40-49 years = 11.5%
- (1 / 3) * log(1 - .115) # 3-year Probability converted to 1-year Rate.
1 - exp(-.04072254 * 1) # 1-year Probability
beta_params(mean = .03990452, sigma = .05)

# 50-59 years = 7.7%
- (1 / 3) * log(1 - .077) # 3-year Probability converted to 1-year Rate.
1 - exp(-.02670868 * 1) # 1-year Probability
beta_params(mean = .02635516,sigma = .05)

# 60-69 years = 5.8%
- (1 / 3) * log(1 - .058) # 3-year Probability converted to 1-year Rate.
1 - exp(-.01991667 * 1) # 1-year Probability
beta_params(mean = .01971964, sigma = .05)

# ===========================================================================================
# Loss to Follow-Up -------------------------------------------------------
# ===========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.
beta_params(mean = .15, sigma = .05)

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