library(rjags)
library(R2jags)
library(bayesplot)
library(parallel)
setwd(dir = "/Users/joshuamusson/Desktop/Analytics/R/Intergrated-CEA-thesis/Evidence_Synthesis/Sub_Models/")
options(mc.cores = detectCores())

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
           rep(p_Age[9], length(30:49)), rep(p_Age[10], length(50:114)))
barplot(p_Age, names.arg = "From Age 15 to 100", ylab = "Probability of HPV+")

N <- 1000
x <- 1:N

plot(density(rpois(n = 1000, lambda = 33.3)))

model_string <- "
model {
 for (i in 1:100) {
  p_Age[i] ~ dgamma(shape, rate)
 }
 
 shape <- pow(m,2) / pow(sd,2)
 rate <-     m / pow(sd,2)
    m ~ dunif(0,100)
    sd ~ dunif(0,100)
}
"
writeLines(text = model_string, con = "AgeInc.txt")

data_JAGS <- list(p_Age = p_Age)

params <- c("shape", "rate")

mod_JAGS <- jags(data = data_JAGS, parameters.to.save = params, model.file = "AgeInc.txt",
                 n.iter = 10000, n.chains = 2, n.burnin = 1000, n.thin = 5)
mod_JAGS
attach.jags(mod_JAGS)

barplot(p_Age, axes = FALSE)
par(new = TRUE)
plot(density(rgamma(n = 10000, shape = mean(shape), rate = mean(rate))), 
     col = "red", lwd = 3, xlab = "From Age 15 to 100", ylab = "Probability of HPV+")

model_String <- model {
  for (i in 1:100) {
    y[i] ~ dgamma
  }
}

weibSurv <- function(t, shape, scale) {
 pweibull(t, shape = shape, scale = scale, lower.tail = F)
}

weibSurv(t = 100, shape = .775, scale = 1 / 40.589)

curve(weibSurv(x, shape = .775, scale = 1 / 0.02463722), from = 0, to = 100,
      ylab = "Survival probability", xlab = "Time", col = "darkblue")
      

40.589 ^ -1
