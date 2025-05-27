## name
dnCDBDP
## title
Character Dependent Birth Death Process
## description
This function simulates a tree under a character-dependent birth-death process.
## details
This function is a flexible simulator that can be used for several phylogenetic models. Specifically,
dnCDBDP allows for character-dependent birrth-death processes meaning that the speciation/extinction rates
depend on the character state of the taxon itself.
 Examples of such models are outlined below:

Multiple State-dependent Speciation Extinction (MuSSE)
This model uses a state-dependent birth-death process to simulate a tree with only anagentic state changes.
Using dnCDBDP to implement MuSSE, a vector of speciation rates for each state can be passed to the lambda argument.
Due to being a vector of rates, dnCDBDP will anly allow anagenetic state changes (along branches), with the
length of the vector corresponding to the number of states.

Cladogenetic State-dependent Speciation Extinction (ClaSSE)
This model allows for cladogenetic state changes (at nodes, ie. state change inducing speciation)  during the birth-death process. 
To implement this, a cladogenetic event map must be passed to the lamda argument which will be a matrix specifiying rates of state changes 
at nodes which will be a seperate matrix for state changes along branches which is specified in the Q argument. See example for implementation.
## authors
Sebastian Hoehna
## see_also
dnCBDSP
fnCladogeneticSpeciationRateMatrix
dnCDCladoBDP
## example
# set up for a two-state ClaSSE model
# set basic starting parameters
root_age ~ dnUniform(0, 2)
rho <- Probability(1/2)

# specifying extinction probabilities for each character
mu_vec[1] := .1
mu_vec[2] := .1

# set up cladogenetic events, probabilities, and number of states
# each element in clado_events describes the cladogenetic event
# for example [0, 0, 1] describes a parent with state 0 and the state of 
# their children, 0 and 1.
clado_events = [[0, 0, 1], [0, 1, 0], [1, 0, 1], [1, 1, 0]]

# set probabilities of each cladogenetic event described above
clado_prob[1] := 1/4
clado_prob[2] := 1/4
clado_prob[3] := 1/4
clado_prob[4] := 1/4

# set total number of states (in this case 0 and 1 are the only states)
num_states = 2

# create cladogenetic rate matrix
clado_matrix = fnCladogeneticProbabilityMatrix(clado_events, clado_prob, num_states)

# set up Q-matrix. rates of state changes along branches 
q_matrix <- matrix([[0, .2], 
                    [.2, 0]])

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
- citation: Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary character's effect on speciation and extinction. Systematic biology, 56(5), 701-710.
  doi: https://doi.org/10.1080/10635150701607033
  url: https://academic.oup.com/sysbio/article/56/5/701/1694265
- citation: FitzJohn, R. G. (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution, 3(6), 1084-1092.
  doi: https://doi.org/10.1111/j.2041-210X.2012.00234.x
  url: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2012.00234.x
- citation: Goldberg, E. E., & Igić, B. (2012). Tempo and mode in plant breeding system evolution. Evolution, 66(12), 3701-3709.
  doi: https://doi.org/10.1111/j.1558-5646.2012.01730.x
  url: https://academic.oup.com/evolut/article/66/12/3701/6851227
