## name
mvMirrorMultiplier
## title
## description
The adaptive mirror multiplier (normal) proposal of Thawornwattana et al. 2017, uses MCMC samples to find posterior mean and variance on the log-scale. After user-defined waiting time, proposes moves (on the log-scale) on opposite side of posterior mean from current location using a normal distribution with the learned posterior standard deviation (scaled by lambda). Before this time, the move uses mu0 as the mean, and lambda as the standard deviation. WARNING: Disabling tuning disables both tuning of proposal variance and learning of empirical mean and variance. To learn the empirical mean and variance without tuning sigma, set adaptOnly=true.
## details
## authors
## see_also
## example
## references
