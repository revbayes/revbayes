## name
dnSBBDP
## title
Sampled speciation birth-death process
## description
Simulates a tree under a birth-death process with a specified sampling rate.
## details
dnSBBDP accepts several arguments to simulate a birth-death process:
- rootAge : Start time for the birth-death process. Accepts a real positive number
- lambda : Vector of speciation rates. Accepts a real positive number
- mu : Vector of extinction rates. Accepts a real positive number
- rho : Taxon sampling probability.  Accepts probability
- taxa : The taxa used for simulation. Accepts Taxon[]
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
                rho           = rho,
                taxa          = TBD)
## references
Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary character's effect on speciation and extinction. Systematic biology, 56(5), 701-710.
