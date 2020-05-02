library(rjags)
library(R2jags)
library(bayesplot)
library(ggplot2)

# ==========================================================================================
# Data on HPV quadrivalent vaccine efficacy ----------------------------------------------
# ==========================================================================================
# The following data will inform sub-model 1: the evidence synthesis for vaccine efficacy.
# All data extracted only includes estimates from the per-protocol pop.

# Function to convert rate to probability assuming constant time:
rateConversion <- function(r, timelength) {
 for (i in 1:timelength) {
  P <- 1 - exp(-r* 1:i)
 }
 print(P)
}

# Study 1: Munoz, et al. 2009 --------------------------------------------------------
# Safety, immunogenicity, and efficacy of quadrivalent human papillomavirus 
# (types 6, 11, 16, 18) recombinant vaccine in women aged 24–45 years: a randomised, 
# double-blind trial. Comparison of vacccine (t = 1) and placebo groups (t = 0).

# RATE OF INFECTION FOR WOMEN t = 0:
# 24-34 = 1.3
# 35-45 = 1.8
# OVERALL = 1.5

# RATE OF INFECTION FOR WOMEN t = 1:
# 24-34 = 0.2
# 35-45 = 0.1
# OVERALL = 0.1

# Study 2: Villa L., et al. 2005 --------------------------------------------------------
# Prophylactic quadrivalent human papillomavirus (types 6, 11, 16, and 18) L1 virus-like 
# particle vaccine in young women: a randomised double-blind placebo-controlled multicentre 
# phase II efficacy trial.
# RATE OF INFECTION FOR WOMEN:
# OVERALL = 0.7
# Total number of observations in vaccine group:

# Study 3: FUTURE II Study Group, 2007 ------------------------------------
# Quadrivalent Vaccine against Human Papillomavirus to Prevent High-Grade Cervical Lesions.

# RATE OF INFECTION FOR WOMEN t = 0:
# OVERALL = .1

# RATE OF INFECTION FOR WOMEN t = 0:
# OVERALL = .3

# ==========================================================================================
# Data on HPV prevalence ----------------------------------------------
# ==========================================================================================
# The following data will inform sub-model 2: the evidence synthesis for HPV prevalence.

 # Proposed model to inform average time to infection with additive age effect:
#  logit(p[i]) <- B0(1 - exp(-r*t)) + p_Age + HIV*I(Yes = 1, No = 0)

# where Indicator shows positive or negative.

# Study 1: Richter K., et al. 2013.  --------------------------------------
# Age-specific prevalence of cervical human papillomavirus infection and cytological 
# abnormalities in women in Gauteng Province. South African Medical Journal, doi-10.7196 
# SAMJ.6514.

# Study population: Women attending public sector primary healthcare clinics for routine 
# gynaecological and non-gynaecological primary healthcare-related reasons.

# Age group by row:
# row 1 = <25
# row 2 = 25-29
# row 3 = 30-34
# row 4 = 35-39 
# row 5 = 40-44
# row 6 = 45-49
# row 7 = 50-54
# row 8 = ≥55
# row 9 = Total

# Total no. of observations:
N <- c(100, 186, 233, 215, 206, 180, 146, 179, 1445)

# HPV+ prevalence:
p_HPV <- c(85, 79.6, 79.4, 72.6, 77.2, 67.8, 64.4, 72.1, 74.6)

# Bind data:
cbind(N, p_HPV)
















