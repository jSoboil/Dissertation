## Copyright Gianluca Baio 2012
lognPar <- function(mu, sdev) {
	sdev.sq <- sdev ^ 2
	mu.log <- log(mu) - .5 * log(1 + sdev.sq / mu ^ 2)
	sdev.sq.log <- log(1 + (sdev.sq / mu ^ 2))
	sigma.log <- sqrt(sdev.sq.log)
	list(mu.log = mu.log, sigma.log= sigma.log)
}

# r_mu.log <- lognPar(mu = r_Age, sdev = r_Age)$mu.log
# r_mu.log
# r_sigma.log <- 1 / lognPar(mu = r_Age, sdev = r_Age)$sigma.log ^ 2
# r_sigma.log

# ==========================================================================================
# Age-specific HPV incidence ----------------------------------------------------
# ==========================================================================================
# Study 1: Sinanovic, E., et al. 2009 -------------------------------------
# The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical
# cancer screening programme in South Africa

# Age-specific incidence:
age_group <- c("15-16", "17", "18", "19", "20", "21", 
                       "22-23", "24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", 
               "55≥")
length(age_group)
group_Age <- c(rep(age_group[1], length(15:16)), rep(age_group[2], 1), 
           rep(age_group[3], 1), rep(age_group[4], 1), rep(age_group[5], 1), 
           rep(age_group[6], 1), rep(age_group[7], length(22:23)), rep(age_group[8], 1),
           rep(age_group[9], length(24:29)), rep(age_group[10], length(30:34)), 
           rep(age_group[11], length(35:39)), rep(age_group[12], length(40:44)), 
           rep(age_group[13], length(45:49)), rep(age_group[14], length(50:54)), 
           rep(age_group[15], length(55:100)))

incidence <- as.numeric(c(.1, .12, .15, .17, .15, .12, rep(.1, 2), .05, rep(.01, 4), rep(.005, 2)))

cbind(incidence, age_group)

r_Age <- c(rep(incidence[1], length(15:16)), rep(incidence[2], 1), 
           rep(incidence[3], 1), rep(incidence[4], 1), rep(incidence[5], 1), 
           rep(incidence[6], 1), rep(incidence[7], length(22:23)), 
           rep(incidence[8], length(24:29)), rep(incidence[9], length(30:34)), 
           rep(incidence[10], length(35:39)), rep(incidence[11], length(40:44)), 
           rep(incidence[12], length(45:49)), rep(incidence[13], length(50:54)), 
           rep(incidence[14], length(55:100)))
r_Age

barplot(r_Age, names.arg = "From Age 15 to 100", 
        main = "Rate of HPV+", sub =  "At constant t-cycle of 1-year")

p_Age <- 1 - exp(-r_Age * 1)

age_group

N <- c(rep(NA, 8), 189, 235, 220, 214, 182, 149, NA)
cbind(age_group, N)
r.Prev <- cbind(rep(NA, length(age_group)), age_group)
r.Prev
n.Prev <- rep(NA, length(p_Age))













