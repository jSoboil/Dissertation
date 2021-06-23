# =====================================================================================
# Gibbs Sampler -----------------------------------------------------------
# =====================================================================================
# Implementing the Gibbs sampling algorithm using a chicken-egg story with unknown hatching 
# probability and invisible unhatched eggs. The first step is to decide on our observed value of X, 
# as well as the constants \lambda, a, and b.
x <- 7
lambda <- 10
a <- 1
b <- 1
# Next state how many iterations to run, and allocate space for the results, creating two vectors 
# p and N of length 104 which we will fill with the simulated draws:
n.iter <- 10^4
p <- rep(0, n.iter)
N <- rep(0, n.iter)

# Initialize p and N to the values 0.5 and 2x, respectively, and then we run the algorithm:
p[1] <- 0.5
N[1] <- 2 * x
for (i in 2:n.iter) {
 p[i] <- rbeta(1, x + a, N[i - 1] - x + b)
 N[i] <- x + rpois(1, lambda * (1 - p[i - 1]))
}

# Add warm-up samples:
p <- p[-(1:(n.iter / 2))]
N <- N[-(1:(n.iter / 2))]

hist(p)
hist(N)

# End file ----------------------------------------------------------------