## name
dnSBBDP
## title
Serially sampled birth-death process
## description
Simulates a tree under a birth-death process with a specified sampling rate.
## details
dnSBBDP simulates a tree under a birth-death process without a character dependent effect.
Additionally, the sampling probability of taxa can be specified using the rho argument of the 
function, allowing for a specfic sampling rate of extinct taxon.

This function is an older version of a birth death process with an explicit sampling rate.
To see a newer version with more funcitonality see: dnBDSTP
## authors
Michael Landis & Sebastian Hoehna
## see_also
dnCDBDP
dnCBDSP
dnBDSTP
## example
# set starting paramters for tree
root_age ~ dnUniform(0, 2)
lambda ~ dnUniform(0, 1)
mu ~ dnUniform(0, 1)
rho := 1/2
# simulate tree
tree ~ dnSBBDP( rootAge       = root_age,
                lambda        = lambda,
                mu            = mu,           
                rho           = rho)
## references
- citation: Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary character's effect on speciation and extinction. Systematic biology, 56(5), 701-710.
  doi: https://doi.org/10.1080/10635150701607033
  url: https://academic.oup.com/sysbio/article/56/5/701/1694265
