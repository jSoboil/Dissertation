# ================================================================================================
# Inverse Transform Sampling ----------------------------------------------
# ================================================================================================
# The following code is used to generate figures that illustrate the inverse transform method:

# Set seed for replication purposes:
set.seed(1000)

# Random samples drawn from uniform distribution, i.e. \theta ~ dunif(0, 1)
x <- runif(1000)
# Plot a histogram of the random draws.
hist(x, col = "white", border = TRUE, lwd = 2, 
     xlab = "Value of the random sample drawn from the uniform distribution", 
     ylab = "Frequency of the randomly drawn value", breaks = 50, 
     main = "Monte Carlo simulation")

# Inverse cdf for exponential distribution:
inv_cdf_dexp <- -log(1 - x) / 3

# Transform samples into an exponential distribution:
hist(inv_cdf_dexp, breaks = 50)

