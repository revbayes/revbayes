## name
dnSBBDP
## title
Serially sampled birth-death process
## description
Simulates a tree under a birth-death process with a specified sampling rate.
## details
dnSBBDP simulates a tree under a birth-death process without a character dependent effect.
Additionally, the sampling probability of taxa can be specified using the rho argument of the 
distribution, allowing for a specfic sampling rate of extinct taxon.

This distribution is an older version of a birth death process with an explicit sampling rate.
To see a newer version with more functionality see: dnBDSTP. dnBDSTP can allow for probability
of death after sampling and additional options for rate interval changes in different rates for the
distribution
## authors
- Michael Landis
- Sebastian Hoehna
## see_also
dnCDBDP
dnCBDSP
dnBDSTP
## example
# set starting paramters for tree
root_age ~ dnUniform(0, 2)
lambda ~ dnUniform(0, 1)
mu ~ dnUniform(0, 1)
rho <- 1
# simulate tree
tree ~ dnSBBDP( rootAge       = root_age,
                lambda        = lambda,
                mu            = mu,           
                rho           = rho)
## references
- citation: Stadler T (2010). Sampling-through-time in birth-death trees. Journal of Theoretical Biology, 267(3):396-404. 
- doi: 10.1016/j.jtbi.2010.09.010
- url: https://www.sciencedirect.com/science/article/pii/S0022519310004765
