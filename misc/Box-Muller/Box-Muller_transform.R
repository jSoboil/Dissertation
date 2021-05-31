library(tidyverse)

# ================================================================================================
# Box-Muller Transform ----------------------------------------------
# ================================================================================================
# The following code is used to generate figures that illustrate the box-muller transform method:

## Random samples drawn from uniform distribution, i.e. \theta ~ dunif(0, 1)
# Set seed for replication purposes:
set.seed(1000)
x1 <- runif(n = 10000)
x2 <- runif(n = 10000)

# Plot a histogram of the random draws.
ggplot() + 
 geom_histogram(aes(x), binwidth = 0.05, colour = "black", fill = "white", alpha = 0.95) + 
 ylab(label = "Frequency") +
 xlab(label = expression(paste(italic(Uniform), " ~ " (0, 1)))) +
 scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
 scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
 theme_light(base_family = "Times New Roman", base_size = 12)

# Transform samples into a normal distribution:
boxMuller_norm <- sqrt(-2 * log(x1)) * cos(2 * pi * x2)

# Plot normal distribution:
ggplot() + 
 geom_density(aes(x = boxMuller_norm), fill = "skyblue", alpha = 0.65) +
 stat_ecdf(aes(boxMuller_norm), colour = "skyblue", alpha = 0.85, lwd = 0.95) + 
 annotate("segment", x = 2.75, xend = 2, y = 0.5, yend = 0.97, colour = "black", 
          size = .4, alpha = 0.75, arrow = arrow()) + 
 annotate("segment", x = 2.75, xend = 1.1, y = 0.5, yend = 0.85, colour = "black", 
          size = .4, alpha = 0.75, arrow = arrow()) +
 annotate("text", x = 2.5, y = .43, label = "CDF of the normal distribution", 
          size = 4.25, family = "Times New Roman") +
 ylab(label = "Probability Density") +
 xlab(label = "Parameter values") +
 theme_linedraw(base_family = "Times New Roman", base_size = 12) + 
 ylim(c(0, 1))
