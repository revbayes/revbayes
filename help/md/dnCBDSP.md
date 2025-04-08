## name
dnCBDSP
## title
Conditional birth-death shift process
## description
Simulates a tree under a birth-death process with shift in birth and death rates.
## details
This function is a flexible simulator that can be used for several phylogenetic models.
dnCBDSP accepts several arguments to simulate a birth-death process:
- rootAge : Start time for the birth-death process. Accepts a real positive number
- rootLambda : Speciation rate at the root the root of tree. Accepts a real positive number
- rootMu : Extinction rate at the root of tree. Accepts a real positive number
- lambda : Prior distribution for speciation rates. Accepts a real positive distribution. Default = NULL
- mu : Prior distribution for extinction rates. Accepts a real positive distribution. Default = NULL
- delta : The rate factor of jumping between rate categories. Accepts real positive numbers
- rho : Taxon sampling probability
- condition : Condition of birth death process. Accepts string. Default = survival. Options: time|survival
- taxa : The taxon names used for initialization. Accepts Taxon[]
## authors
## see_also
dnCDBDP
dnSBBDP
## example
# set distributions for tree
root_age ~ dnUniform(0, 2)
root_lambda ~ dnUniform(0, 1)
root_mu ~ dnUniform(0, 1)
sampling_prob := 1/2
# simulate tree
tree ~ dnCBDSP( rootAge           = root_age,
                rootLambda        = root_lambda,
                rootMu            = root_mu,
                delta             = .2,
                rho               = sampling_prob,
                condition         = "survival"
		taxa              = TBD)
## references
Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary character's effect on speciation and extinction. Systematic biology, 56(5), 701-710.

