library(tidyverse)

# ================================================================================================
# Inverse Transform Sampling ----------------------------------------------
# ================================================================================================
# The following code is used to generate figures that illustrate the inverse transform method:

## Random samples drawn from uniform distribution, i.e. \theta ~ dunif(0, 1)
# Set seed for replication purposes:
set.seed(1000)
x <- runif(n = 1000)

# Plot a histogram of the random draws.
ggplot() + 
 geom_histogram(aes(x), binwidth = 0.05, colour = "black", fill = "white", alpha = 0.95) + 
 ylab(label = "Frequency") +
 xlab(label = expression(paste(italic(Uniform), " ~ " (0, 1)))) +
 scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
 scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
 theme_light(base_family = "Times New Roman", base_size = 12)

# Transform samples into a exponential distribution:
inv_cdf_dexp <- qexp(x)
# Transform samples into a normal distribution:
inv_cdf_norm <- qnorm(p = x)

# Plot of both exponential and normal distribution:
ggplot() + 
 geom_density(aes(x = inv_cdf_norm), fill = "blue", alpha = .65) +
 geom_density(aes(x = inv_cdf_dexp), fill = "orange", alpha = .45) + 
 scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
 scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5)) + 
 ylab(label = "Probability Density") +
 xlab(label = "Parameter values") +
 theme_minimal(base_family = "Times New Roman", base_size = 12)

# Just to show how this is equivalent to using the rnorm function built into R:
x <- runif(n = 10000)
x <- qnorm(p = x, mean = 0.8, sd = 1)
plot(density(x))
plot(ecdf(x))
# ... now using the random distribution function in R...
y <- rnorm(n = 10000, mean = 0.8, sd = 1)
plot(ecdf(y))

