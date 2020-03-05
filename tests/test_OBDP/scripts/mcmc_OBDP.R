################################################################################
#
# RevBayes Validation Test: Fossilized birth-death process
#
# Model: Tree is drawn from a constant-rate fossilized birth-death process.
#
#
# authors: Walker Pett
#
################################################################################

seed(12345)

#######################
# Reading in the Data #
#######################

### Read in taxon data


# set my move index
mvi = 0
taxa <- readTaxonData("data/bears_taxa.tsv")
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)
rho <- 1.0
rm <- 0.0
cond <- "time"
omega <- 0.0
origin_time ~ dnUnif(37.0, 55.0)




 # create some moves that change the stochastic variables
 # all moves are sliding proposals but you could use scaling proposals for the rates too
 moves[mvi++] = mvScale(lambda,lambda=1,weight=1)
 moves[mvi++] = mvScale(lambda,lambda=0.1,weight=1)
 moves[mvi++] = mvScale(lambda,lambda=0.01,weight=1)

 moves[mvi++] = mvScale(mu,lambda=1,weight=1)
 moves[mvi++] = mvScale(mu,lambda=0.1,weight=1)
 moves[mvi++] = mvScale(mu,lambda=0.01,weight=1)

 moves[mvi++] = mvScale(psi,lambda=1,weight=1)
moves[mvi++] = mvScale(psi,lambda=0.1,weight=1)
moves[mvi++] = mvScale(psi,lambda=0.01,weight=1)


### Define the tree-prior distribution as the fossilized birth-death process ###
obd_dist = dnOBDP(originAge=origin_time, lambda=lambda, mu=mu, psi=psi,omega=omega, dn_time_points=v(0.0), rho=rho, removalPr=rm, cond=cond, taxa=taxa)
print("OK dn")

# The will be a random variable of a constrained topology distribution that is governed by the obd #
# this distribution will generate obd trees that match the monophyly constraints defined above #
clade_ursinae = clade("Melursus_ursinus", "Ursus_arctos", "Ursus_maritimus",
                      "Helarctos_malayanus", "Ursus_americanus", "Ursus_thibetanus",
                      "Ursus_abstrusus", "Ursus_spelaeus")
constraints = v(clade_ursinae)


obd_tree ~ dnConstrainedTopology(obd_dist, constraints=constraints)
print("OK tree")

moves[mvi++] = mvFNPR(obd_tree, weight=1.0)
moves[mvi++] = mvTipTimeSlideUniform(obd_tree, weight=1.0)
moves[mvi++] = mvNodeTimeSlideUniform(obd_tree, weight=1.0)
moves[mvi++] = mvRootTimeSlideUniform(obd_tree, origin_time, weight=1.0)


#############
# THE Model #
#############

# We define our model.
# We can use any node of our model as a handle, here we chose to use the rate matrix.
mymodel = model(obd_tree)


monitors[1] = mnStochasticVariable(filename="output/mcmc_OBDP.out",printgen=100)
monitors[2] = mnFile(filename="output/mcmc_OBDP.trees",printgen=100,obd_tree)

mymcmc = mcmc(mymodel, monitors, moves, moveschedule="single")


mymcmc.run(generations=100000)


# Read in the tree trace and construct the maximum clade credibility (MCC) tree #
trace = readTreeTrace("output/mcmc_OBDP.trees")

# Summarize tree trace and save MCC tree to file
mccTree(trace, file="output/mcmc_OBDP_mcc.tre" )


you may want to quit RevBayes now
q()
