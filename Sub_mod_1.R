library(rjags)
library(R2jags)
library(bayesplot)
library(ggplot2)

# r <- c(rep(85, 7), 149, 185, 157, 159, 122, 94, 129)
r <- c(rep(85, 7))
# n <- c(rep(100, 7), 186, 233, 215, 206, 180, 146, 179)
n <- c(rep(100, 7))
# p_Age <- c(0.1, 0.12, 0.15, 0.17, 0.15, 0.12, 0.1, 0.1, 0.1, 
#           0.1, 0.05, 0.01, 0.005, 0.005)
p_Age <- c(0.1, 0.12, 0.15, 0.17, 0.15, 0.12, 0.1)

model_String <- "
 # Sub-model 1: Population level incidence for Ages 16-23:
   model {
    for (i in 1:length(r)) {
     # Binomial Likelihood on 
     # postive outcome:
     r[i] ~ dbin(p[i], n[i])
  
     # Link function between
     # Age and Outcome:
     p[i] <- p_Age[i]
   
     # Normal Likelihood on 
     # Age:
     p_Age[i] ~ dnorm(mu, prec)
   }
  
 # Priors on mu and sd:
 mu ~ dnorm(0, .1)
 sd ~ dunif(0, 10)
 tau.sq <- sd * sd
 prec <- 1 / tau.sq
  
 # Priors on LSIL and
 # HSIL
 LSIL ~ dbeta(.2, .8)
 HSIL ~ dbeta(.35, .9)
 
}
"
writeLines(text = model_String, con = "test.txt")


data_JAGS <- list(r = r, n = n, p_Age = p_Age)
data_JAGS

inits <- list(
 list(LSIL = c(.1), HSIL = c(.1)),
 list(LSIL = c(.3), HSIL = c(.4)) 
 )

params <- c("LSIL", "HSIL", "mu", "tau.sq")

jags_Model <- jags(data = data_JAGS, parameters.to.save = params, model.file = "test.txt", 
                   inits = inits, n.chains = 2, n.iter = 10000, n.burnin = 1000, n.thin = 18)
jags_Model

posterior <- as.array(jags_Model$BUGSoutput$sims.array)
dimnames(posterior)

color_scheme_set("viridisB")
theme_set(theme_minimal())
mcmc_trace(posterior, pars = params, 
           facet_args = list(ncol = 1, strip.position = "left"))

color_scheme_set("viridisA")
theme_set(theme_minimal())
mcmc_trace(posterior[,, 1:4], window = c(100, 150), size = 1) + 
  panel_bg(fill = "white", color = NA) +
  legend_move("top")

color_scheme_set("mix-teal-pink")
mcmc_dens_overlay(posterior, params)

color_scheme_set("pink")
mcmc_pairs(posterior, pars = params,
           off_diag_args = list(size = 1.5))

color_scheme_set("mix-blue-brightblue")
mcmc_acf(posterior, pars = params, lags = 50)

color_scheme_set("mix-teal-pink")
mcmc_hist(x = posterior, pars = params, binwidth = .001)


