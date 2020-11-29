# ============================================================================================
# Metropolis-Hastings -----------------------------------------------------
# ============================================================================================
# The Metropolis-Hastings algorithm implemented on a Normal-Normal model.

# First, we must choose our observed value of Y and decide on values for the constants \sigma, 
# \mu, \tau.

# Likelihood normal model:
y <- 3
sigma <- 1
# Prior normal:
mu <- 0
tau <- 2

# We also need to choose the standard deviation of the proposals for step 1 of the algorithm, for 
# this problem, we let d = 1. We set the number of iterations to run, and we allocate a vector 
# \theta of length 104 which we will fill with our simulated draws:
d <- 1
n.iter <- 10^4
theta <- rep(0, n.iter)

# Now for the main loop. We initialize ✓ to the observed value y, then run the algorithm:
theta[1] <- y
for (i in 2:n.iter) {
 theta.p <- theta[i - 1] + rnorm(1, 0, d)
 
 r <- (dnorm(y, theta.p, sigma) * dnorm(theta.p, mu, tau)) / 
  (dnorm(y, theta[i - 1], sigma) * dnorm(theta[i - 1], mu, tau))
 
 flip <- rbinom(1, 1, min(r, 1))
 
 theta[i] <- if (flip == 1) theta.p else theta[i - 1]
}
# Warmup:
theta <- theta[ -(1:(n.iter / 2))]

hist(theta)
mean(theta)
#  Note the use of ratios in the MCMC algorithm.

# The proposed value of \theta is theta.p, which equals the previous value of \theta plus a Normal
# random variable with mean 0 and standard deviation d. theta.p is playing the role of x' and 
# theta[i-1] is playing the role of x. The coin flip to determine whether to accept or reject the 
# proposal is flip, which is a coin flip with probability min(r,1) of Heads (encoding Heads as 1 
# and Tails as 0). Finally, we set theta[i] equal to the proposed value if the coin flip lands 
# Heads, and keep it at the previous value otherwise. theta <- theta[-(1:(niter/2))] provides a 
# warmup.

# End file ----------------------------------------------------------------