###############################################################################
#
# RevBayes Validation Test: Occurrence birth-death process
#
# Model: Tree is drawn from a constant-rate fossilized birth-death process with occurrences.
#
#
# authors: Antoine Zwaans (from Walker Pett's FBDP test function)
#
################################################################################

seed(12345)

#######################
# Reading in the Data #
#######################

### Read in taxon data


# set my move index
mvi = 1
taxa <- readTaxonData("data-epi-test/data_taxa.csv")
sequences <- readDiscreteCharacterData("data-epi-test/data_seq.nex")
occurrence_ages <- readMatrix(file="data-epi-test/data_occurrences.csv", delimiter="; ")[1]



rho <- 1.0
rm <- 0.0
cond <- "time"
origin_time ~ dnUnif(7.7, 12.0)
N <- 15
Mt <- TRUE

#Add Missing Taxa
sequences.addMissingTaxa( taxa )

# Create some moves that change the stochastic variables
# All moves are sliding proposals but you could use scaling proposals for the rates too
lambda ~ dnExp(10)
moves[mvi++] = mvScale(lambda,lambda=0.1,tune=TRUE,weight=1.0)
mu ~ dnExp(10)
moves[mvi++] = mvScale(mu,lambda=0.1,tune=TRUE,weight=1.0)
psi ~ dnExp(10)
moves[mvi++] = mvScale(psi,lambda=0.1, tune=TRUE,weight=1.0)
omega ~ dnExp(10)
moves[mvi++] = mvScale(omega,lambda=0.1, tune=TRUE,weight=1.0)

### Define the tree-prior distribution as the fossilized birth-death process ###
obd_tree ~ dnOBDP(originAge=origin_time, lambda=lambda, mu=mu, psi=psi,omega=omega, maxHiddenLin=N, rho=rho, removalPr=rm, cond=cond,occurrence_ages=occurrence_ages,useMt=Mt, taxa=taxa)

moves[mvi++] = mvFNPR(obd_tree, weight=10.0)
moves[mvi++] = mvCollapseExpandFossilBranch(obd_tree, origin_time, weight=6.0)
moves[mvi++] = mvNodeTimeSlideUniform(obd_tree, weight=2.0)
moves[mvi++] = mvRootTimeSlideUniform(obd_tree, origin_time, weight=1.0)


#strict clock model with lambda = 0.005
branch_rates ~ dnExponential(250.0)
moves[mvi++] = mvScale(branch_rates, lambda=0.01, tune=TRUE)


#nucleotide evolution HKY model
pi_prior <- v(1,1,1,1)
pi ~ dnDirichlet(pi_prior)

moves[mvi++] = mvBetaSimplex(pi, weight=2)
moves[mvi++] = mvDirichletSimplex(pi, weight=1)

kappa ~ dnLognormal(0.0, 1.0)
moves[mvi++] =  mvScale(kappa)

Q_epi := fnHKY(kappa,pi)

alpha_epi ~ dnExponential( 1.0 )
moves[mvi++] =  mvScale(alpha_epi, lambda=0.1, tune=TRUE)

rates_epi := fnDiscretizeGamma( alpha_epi, alpha_epi, 4 )

phySeq ~ dnPhyloCTMC(tree=obd_tree, Q=Q_epi, siteRates=rates_epi, branchRates=branch_rates, type="DNA")
phySeq.clamp(sequences)


#############
# THE Model #
#############
mymodel = model(obd_tree)

# We define our model.
# We can use any node of our model as a handle, here we chose to use the rate matrix.
monitors[1] = mnStochasticVariable(filename="out-mini-epi/mcmc_OBDP_epi_new.out", printgen=10)
monitors[2] = mnFile(filename="out-mini-epi/mcmc_OBDP_epi_new.trees", printgen=100,obd_tree)
monitors[3] = mnScreen(printgen=100)
print("my model ok")
#mymcmc = mcmc(mymodel, monitors, moves, moveschedule="single")
mymcmc = mcmc(mymodel, monitors, moves)

mymcmc.run(generations=5000, tuningInterval=100)

# check the performance of the MCMC/moves
mymcmc.operatorSummary()

# Read in the tree trace and construct the maximum clade credibility (MCC) tree #
trace = readTreeTrace("out-mini-epi/mcmc_OBDP_epi_new.trees")

# Summarize tree trace and save MCC tree to file
mccTree(trace, file="out-mini-epi/mcmc_OBDP_epi_new.tre" )


# you may want to quit RevBayes now
q()
