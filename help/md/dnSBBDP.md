## name
dnSBBDP
## title
Sampled-speciation birth-death process
## description
## details
## authors
Michael Landis
Sebastian Hoehna
## see_also
dnCBDSP
dnCDBDP
## example
# set parameters for the process
taxa <- [taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E")]
root_age ~ dnUniform(0, 2)
lambda ~ dnUniform(0, 1)
mu ~ dnUniform(0, 1)
rho <- 1
# simulate tree
tree ~ dnSBBDP(rootAge = root_age,
               lambda  = lambda,
               mu      = mu,           
               rho     = rho,
               taxa    = taxa)
## references
