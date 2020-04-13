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

observed_phylogeny <- readTrees("output_macro_test_taxa/mcmc_OBDP_macro_test.tre")[1]
taxa <- observed_phylogeny.taxa()

# set my move index
mvi = 0


##############
# Tree model #
##############


# specify a prior on the origin age
#origin_time ~ dnUnif(37.0, 55.0)

# first we create the standard deviation of the rates between intervals
# draw the sd from an exponential distribution
speciation_sd ~ dnExponential(1.0)
moves[++mvi] = mvScale(speciation_sd,weight=10.0)

extinction_sd ~ dnExponential(1.0)
moves[++mvi] = mvScale(extinction_sd,weight=10.0)


# create a random variable at the present time
log_speciation[1] ~ dnUniform(-10.0,10.0)
log_extinction[1] ~ dnUniform(-10.0,10.0)


# apply moves on the rates
moves[++mvi] = mvSlide(log_speciation[1], weight=5)
moves[++mvi] = mvSlide(log_extinction[1], weight=5)


speciation[1] := exp( log_speciation[1] )
extinction[1] := exp( log_extinction[1] )

NUM_INTERVALS = 10
for (i in 1:NUM_INTERVALS) {
    index = i+1

    # specify normal priors (= Brownian motion) on the log of the rates
    log_speciation[index] ~ dnNormal( mean=log_speciation[i], sd=speciation_sd )
    log_extinction[index] ~ dnNormal( mean=log_extinction[i], sd=extinction_sd )

    # apply moves on the rates
    moves[++mvi] = mvSlide(log_speciation[index], weight=5)
    moves[++mvi] = mvSlide(log_extinction[index], weight=5)

    # transform the log-rate into actual rates
    speciation[index] := exp( log_speciation[index] )
    extinction[index] := exp( log_extinction[index] )

}

moves[++mvi] = mvVectorSlide(log_speciation, weight=20)
moves[++mvi] = mvVectorSlide(log_extinction, weight=20)

moves[++mvi] = mvShrinkExpand( log_speciation, sd=speciation_sd, weight=20 )
moves[++mvi] = mvShrinkExpand( log_extinction, sd=extinction_sd, weight=20 )


interval_times <- observed_phylogeny.rootAge() * (1:NUM_INTERVALS) / (NUM_INTERVALS) * 0.8
print("speciation rates = ")
print(speciation)
print("interval times = ")
print(interval_times)

origin_time = observed_phylogeny.rootAge()
### Define the tree-prior distribution as the fossilized birth-death process ###
fbd_tree ~ dnOBDP(origin=origin_time, lambda=speciation, mu=extinction, psi=speciation, timeline=interval_times, rho=1.0, taxa=taxa, initialTree=observed_phylogeny,useMt
=FALSE)


# The will be a random variable of a constrained topology distribution that is governed by the FBD #
# this distribution will generate FBD trees that match the monophyly constraints defined above #


moves[mvi++] = mvFNPR(fbd_tree, weight=1.0)
#moves[mvi++] = mvCollapseExpandFossilBranch(fbd_tree, origin_time, weight=1.0)
moves[mvi++] = mvTipTimeSlideUniform(fbd_tree, weight=1.0)
moves[mvi++] = mvNodeTimeSlideUniform(fbd_tree, weight=1.0)
#moves[mvi++] = mvRootTimeSlideUniform(fbd_tree, origin_time, weight=1.0)

num_samp_anc := fbd_tree.numSampledAncestors()
num_tips := fbd_tree.ntips()

#############
# THE Model #
#############

# We define our model.
# We can use any node of our model as a handle, here we chose to use the rate matrix.
mymodel = model(fbd_tree)


monitors[1] = mnStochasticVariable(filename="output_pOBDP/mcmc_pOBDP.out",printgen=1)
monitors[2] = mnFile(filename="output_pOBDP/mcmc_pOBDP.trees",printgen=1,fbd_tree)
monitors[3] = mnScreen(printgen=1, num_samp_anc, num_tips)

mymcmc = mcmc(mymodel, monitors, moves, moveschedule="single")

mymcmc.run(generations=1)

# check the performance of the MCMC/moves
mymcmc.operatorSummary()

# Read in the tree trace and construct the maximum clade credibility (MCC) tree #
#trace = readTreeTrace("output_FBDRP/mcmc_FBDRP.trees")

# Summarize tree trace and save MCC tree to file
#mccTree(trace, file="output_FBDRP/mcmc_FBDRP_mcc.tre")


# you may want to quit RevBayes now
q()
