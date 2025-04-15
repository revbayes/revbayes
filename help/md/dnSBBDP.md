## name
dnSBBDP
## title
Sampled speciation birth-death process
## description
Simulates a tree under a birth-death process with a specified sampling rate.
## details
dnSBBDP simulates a tree under a birth-death process without a character dependent effect.
Additionally, the sampling probability of taxa can be specified using the rho argument of the 
function, allowing for a specfic sampling rate of extinct taxon.

## authors
Michael Landis & Sebastian Hoehna
## see_also
dnCDBDP
dnCBDSP
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
