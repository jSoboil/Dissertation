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

# ==========================================================================================
# Age-specific all cause mortality --------------------------------------------------
# ==========================================================================================
# Import ASSA mortality table:
mort_data <- read_excel("data-raw/mortality tables.xls", 
                        sheet = "ASSA data", range = "A3:H94")
mort_data
# Total Pop. divided by total deaths:
v_p_mort_All <- mort_data[1:86, 3] / mort_data[1:86, 2]
v_p_mort_All

# Probability of mortality Cervical Cancer from HPV:
v_p_mort_HPV <- mort_data[1:86, 8]
colnames(v_p_mort_HPV) <- "Death_HPV"
v_p_mort_HPV

# Probability of mortality less HPV:
v_p_mort_lessHPV <- v_p_mort_All - v_p_mort_HPV
v_p_mort_lessHPV <- cbind(mort_data$Age[1:86], v_p_mort_lessHPV)
colnames(v_p_mort_lessHPV) <- c("Age", "Death_less.HPV")
v_p_mort_lessHPV

# ===========================================================================================
# Collated studies on vaccine efficacy -------------------------------------------------
# ===========================================================================================

# Study 1: Apter et al. (2015)   ------------------------------------------
# Efficacy of Human Papillomavirus 16 and 18 (HPV-16/18) AS04- Adjuvanted Vaccine against 
# Cervical Infection and Pre-cancer in Young Women: Final Event-Driven Analysis of the 
# Randomized, Double- Blind PATRICIA Trial.

# Continent: 

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
# p = 1 - exp(-r * t)
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


# ==========================================================================================
# From Infection to LSIL -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

1 - exp(-.2 * 3)

# ==========================================================================================
# From Infection to HSIL -----------------------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.


# ==========================================================================================
# From LSIL to Infection or Normal ----------------------------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Ages 15-34:
1 - exp(-.65 * 6)
# Proportion reverting to Normal is .9:
1 - ((0.9797581 - 0.9797581 * .9) + 0.9797581 * .9)

# Ages ≥ 35:
1 - exp(-.4 * 6)
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

# Probability of progression every 10 years = .4
# Convert to yearly rate
- (1 / 10) * log(1 - .4)
# Convert to annual progression probability to Stage I Cancer:
1 - exp(-0.05108256 * 1)

# ==========================================================================================
# From Stage I Cancer to Treatment or Stage II Cancer ---------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Annual probability from Stage I to Stage II:
(1 - exp(-0.9 * 4))

# Annual probability of symptoms from Stage I to Stage I Treatment:
.15

# Yearly probability from Stage I to Stage I:
# 1 - (StageI_to_StageII + StageI_to_Death + StageI_to_Treatment)

# To check probability laws, uncomment the line below:
# 0.710708 + (1 - exp(-.15 * 1) + .15)

# ==========================================================================================
# From Stage II Cancer to Treatment or Stage III Cancer -------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Probability of progression every 3 years = .9
# Convert to yearly rate

# Yearly probability from Stage I to Stage II:
1 - exp(-0.9 * 3)

# Annual probability of symptoms from Stage II to Treatment:
0.225

# Yearly probability from Stage II to Stage II:
# 1 - (StageII_to_StageIII + StageII_to_Death + StageII_to_Treatment)

# ==========================================================================================
# From Stage III Cancer to Treatment or Stage IV Cancer -------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

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
# 1 - (StageIII_to_StageIV + StageIII_to_Death + StageIII_to_Treatment)

# ==========================================================================================
# From Stage IV Cancer to Treatment -------------------
# ==========================================================================================
# Informative prior -------------------------------------------------------
# 1. Sinanovic, E., et al. 2009. The potential cost-effectiveness of adding a human 
# papillomavirus vaccine to the cervical cancer screening programme in South Africa.

# Yearly probability of symptoms from Stage IV to Treatment:
.9

# Yearly probability from Stage IV to Stage IV:
# 1 - (StageIV_to_Death + StageIV_to_Treatment)

# ==========================================================================

# ===========================================================================================
# Stage I 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.9688
# Die = 1 - 0.9688

# Year 2
# Survive = 0.9525
# Die = 1 - 0.9525

# Year 3
# Survive = 0.9544
# Die = 1 - 0.9544

# Year 4
# Survive = 0.9760
# Die = 1 - 0.9760

# Year 5
# Survive = 0.9761
# Die = 1 - 0.9761

# 5 year survival:
((0.9688) * (0.9525) * (0.9544) * (0.9760) * (0.9761))

# ===========================================================================================
# Stage II 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.9066
# Die = 1 - 0.9066

# Year 2
# Survive = 0.8760
# Die = 1 - 0.8760

# Year 3
# Survive = 0.9225
# Die = 1 - 0.9225

# Year 4
# Survive = 0.9332
# Die = 1 - 0.9332

# Year 5
# Survive = 0.9604
# Die = 1 - 0.9604

# 5 year survival:
((0.9066) * (0.8760) * (0.9225) * (0.9332) * (0.9604))

# ===========================================================================================
# Stage III 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.7064
# Die = 1 - 0.7064

# Year 2
# Survive = 0.7378
# Die = 1 - 0.7378

# Year 3
# Survive = 0.8610
# Die = 1 - 0.8610

# Year 4
# Survive = 0.9231
# Die = 1 - 0.9231

# Year 5
# Survive = 0.9142
# Die = 1 - 0.9142

# 5 year survival:
((0.7064) * (0.7378) * (0.8610) * (0.9231) * (0.9142))

# ===========================================================================================
# Stage IV 5-year Survival -------------------------------------------------
# ===========================================================================================
# Year 1:
# Survive = 0.3986
# Die = 1 - 0.3986

# Year 2
# Survive = 0.4982
# Die = 1 - 0.4982

# Year 3
# Survive = 0.7638
# Die = 1 - 0.7638

# Year 4
# Survive = 0.8652
# Die = 1 - 0.8652

# Year 5
# Survive = 0.8592
# Die = 1 - 0.8592

# 5 year survival:
((0.3986) * (0.4982) * (0.7638) * (0.8652) * (0.8592))

# ==========================================================================================
# Risk Increase (HIV+) ------------------------------------------------------
# ==========================================================================================
# Note: although desireble, we did not add a risk parameter to the final model and chose to 
# follow the core assumptions of the original model.

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