## name
dnCDBDP
## title
Character Dependent Birth Death Process
## description
This function simulates a tree under a character-dependent birth-death process.
## details
This function is a flexible simulator that can be used for several phylogenetic models.
dnCDBDP accepts several arguments to simulate a birth-death process:
- rootAge : Start time for the birth-death process. Accepts a Real positive number
- speciationRates/lamda/cladoEventMap : Vector of speciation rate if anagenetic-only model or cladogenetic event map.
- extinctionRates/mu : Vector of extinction rates. Accepts real positive numbers
- psi/phi : Vector of serial sample rates. Accepts real positive numbers. Default = NULL
- Q : The rate matrix of jumping between categories. Default = NULL
- delta : The rate factor of jumping between categories. Accepts real positive numbers. Default = 1
- pi : Root state frequencies. Accepts simplex. Default = NULL
- rho : Taxon sampling probability. Default = 1
- condition : Condition of birth death process. Accepts string. Default = time. Options: time|survival
- nTimeSlices : Number of time slices for numeric ODE. Accepts real positive number. Default = 500
- simulateCondition : Conditions under which to simulate. Accepts string. Default = startTime. Options: startTime|numTips|tipStates|tree
- minNumLineages : Minimum number of lineages to simulate; applied under startTime condition. Accepts a natural number. Default = 0
- maxNumLineages : Maximum number of lineages to simulate; applied under startTime condition. Accepts a natural number. Default = 500
- exactNumLineages : Exact number of lineages to simulate; applied under numTips and tipStates conditions. Accepts a natural number. Default = 100
- maxTime : Maximum time for lineages to coalesce when simulating; applied under the numTips and tipStates condition. Accepts a real positive number. Default = 1000.
- pruneExtinctLineages : Should simulation prune extinct lineages? Accepts boolean. Default = TRUE.
- allowRateShiftsAtExtinctLineages : Should we allow rate shifts to occur on extinct lineages?. Accepts boolean. Default = TRUE.
## authors
Sebastian Hoehna
## see_also
dnCBDSP
fnCladogeneticSpeciationRateMatrix
dnCDCladoBDP
## example
# set basic starting parameters
root_age ~ dnUniform(0, 2)
rho := Probability(1/2)

# set up cladogenetic events, probabilities, and number of states
clado_events = [[0, 0, 1], [0, 1, 0], [1, 0, 1], [1, 1, 0]]
clado_prob[1] := 1/4
clado_prob[2] := 1/4
clado_prob[3] := 1/4
clado_prob[4] := 1/4
num_states = 2

# create cladogenetic rate matrix
clado_matrix = fnCladogeneticProbabilityMatrix(clado_events, clado_prob, num_states)

# specifiying extinction probabilites
mu_vec[1] := .1
mu_vec[2] := .1

# set up Q-matrix
q_matrix <- matrix([[0, .2], [.2, 0]])

# set pi, root state freq
pi <- simplex([1, 2])

# basic use of the function
timetree ~ dnCDBDP( rootAge           = root_age,
                    speciationRates   = clado_matrix,
                    extinctionRates   = mu_vec,
                    Q                 = q_matrix,
                    pi                = pi,
                    rho               = rho,
                    condition         = "time")

## references
Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary character's effect on speciation and extinction. Systematic biology, 56(5), 701-710.

FitzJohn, R. G. (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution, 3(6), 1084-1092.

Goldberg, E. E., & Igić, B. (2012). Tempo and mode in plant breeding system evolution. Evolution, 66(12), 3701-3709.
