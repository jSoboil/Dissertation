library(rjags)
library(R2jags)
library(bayesplot)
library(ggplot2)

r <- c(rep(85, 7), 149, 185, 157, 159, 122, 94, 129)

n <- c(rep(100, 7), 186, 233, 215, 206, 180, 146, 179)

p_Age <- c(0.1, 0.12, 0.15, 0.17, 0.15, 0.12, 0.1, 0.1, 0.1, 
           0.1, 0.05, 0.01, 0.005, 0.005)

# Population prevalence:
# ... based on eq. from (Fleurence and Hollenbeak, 2007).
# From ICO estimated 3.2 of population have HPV at any given time...
pop_Prev <- - 1 * log(1 - .032)
pop_Prev

model_String <- " # Sub-model 1: Population prevalence given data:
   model {
    for (i in 1:length(r)) {
     # Binomial Likelihood on 
     # postive outcome:
     r[i] ~ dbin(pop_Prev[i], n[i])
     
     pop_Prev[i] ~ dnorm(.032, prec)
     
    }
 # Prior on prec:
 tau ~ dunif(0, 100)
 tau.sq <- tau * tau
 prec <- 1 / tau.sq
 # Priors on LSIL and
 # HSIL
 LSIL ~ dbeta(.2, .8)
 HSIL ~ dbeta(.35, .9)
 
}
"
writeLines(text = model_String, con = "test.txt")


data_JAGS <- list(n = n, r = r)
data_JAGS

inits <- list(
 list(LSIL = c(.1), HSIL = c(.1)),
 list(LSIL = c(.3), HSIL = c(.4)) 
 )

params <- c("LSIL", "HSIL", "pop_Prev")

jags_Model <- jags(data = data_JAGS, parameters.to.save = params, model.file = "test.txt", 
                   inits = inits, n.chains = 2, n.iter = 10000, n.burnin = 1000, 
                   n.thin = 18)
jags_Model

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = c("HSIL", "LSIL", "pop_Prev"), 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[,, 1:4], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, pars = c("HSIL", "LSIL"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("HSIL", "LSIL"),
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = c("HSIL", "LSIL"), lags = 50)

color_scheme_set("mix-teal-pink")
mcmc_hist(x = posterior, pars = c("HSIL", "LSIL"), 
          binwidth = .035)


