###############################################################################
#
# RevBayes Validation Test: Occurrence birth-death process
#
# Model: Tree is drawn from a constant-rate fossilized birth-death process with occurrences.
#
#
# authors: Antoine Zwaans + Jérémy Andréoletti
#
################################################################################
seed(12345)
#######################
# Reading in the Data #
#######################

### Read in taxon data
tree <- readTrees("data_OBDP2_simple/tree.tre")[1]
taxa <- tree.taxa()
occurrence_ages <- readMatrix(file="data-OBDP2-simple/occurrences.csv", delimiter=";")[1]

#set rate values
speciation <- v(1.0,1.0,1.0)
extinction <- 0.9
rho <- 0.2
sampling <- 0.25
omega <- 0.5
rm <- 0.3
N <- 50
cond <- "survival2"
origin_time <- 9
Mt <- TRUE

#set timeline
timeline <- v(1,2)

print(speciation)
### Define the tree-prior distribution as the piecewise constant birth-death process ###
obd_tree ~  dnOBDP2(originAge=origin_time,
                    timeline=timeline,
                    lambda=speciation,
                    mu=extinction,
                    omega=omega,
                    phi=sampling,
                    Phi=rho,
                    cond=cond,
                    r=rm,
                    occurrence_ages=occurrence_ages,
                    taxa=taxa,
                    useMt=Mt,
                    verbose=TRUE,
                    initialTree=tree)


moves[1] = mvFNPR(obd_tree, weight=1.0)
#############
# THE Model #
#############
# We define our model.
# We can use any node of our model as a handle, here we chose to use the rate matrix.
mymodel = model(obd_tree)
monitors[1] = mnScreen(printgen=1)
print("my model ok")
mymcmc = mcmc(mymodel, monitors, moves, moveschedule="single")
mymcmc.run(generations=1, tuningInterval=1)
# you may want to quit RevBayes now
q()
