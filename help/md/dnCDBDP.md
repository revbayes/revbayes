## name
dnCDBDP
## title
Character-dependent birth-death process
## description
Simulates a tree under a multi-type birth-death process.
## details
This distribution is a flexible simulator that can be used for several
lineage-heterogeneous diversification models. Specifically, `dnCDBDP` allows
the tree to be divided into multiple segments (also referred to as "regimes"
or "types"; Barido-Sottani et al. 2020) separated by discrete shifts in birth
and death rates, such that each segment draws its rate vector from a finite
number of categories. These categories can (but need not) correspond to the
states of an observed discrete character, in which case the shifts correspond
to inferred state transitions: hence the characterization of `dnCDBDP` as
a "character-dependent" model.

Applications of this multi-type model include:

1. Multiple State-dependent Speciation Extinction (MuSSE) (Maddison et al.
   2007; FitzJohn 2012)

This model uses a state-dependent birth-death process to simulate a tree with
only anagenetic state changes, i.e., changes taking place along the branches
of the tree. When `dnCDBDP` is used to implement MuSSE, the `lambda` argument
represents a vector of speciation rates for each state, and its length is
therefore equal to the number of states. The rates at which the character
transitions from one state to another are specified by the matrix `Q` (with as
many rows and columns as there are states), and the extinction rates for each
state are specified using the vector `mu` (with as many elements as there are
states). A draw from the distribution needs to be fixed both to a previously
inferred tree (using `.clamp()`), and to a matrix recording which state is
observed at a given tip (using `.clampCharData()`).

2. Cladogenetic State-dependent Speciation Extinction (ClaSSE) (Goldberg &
   Igić 2012)

This model additionally allows for cladogenetic state changes, i.e., changes
that take place at nodes, corresponding to the assumption that state change
either induces or immediately follows speciation. When `dnCDBDP` is used to
implement ClaSSE, the `lambda` argument corresponds to a cladogenetic event
map. A cladogenetic event map is a matrix specifying the speciation rates
associated with the different character state triplets that can be observed
for the parent and its two children at a given node (= cladogenetic event).
As before, anagenetic state change is still allowed, with its rates described
by the `Q` matrix, and each state is associated with a distinct extinction
rate specified in the `mu` vector. A draw from the distribution again has to 
be clamped both to a tree and to a character matrix. See the example below
for implementation.

3. Branch-specific Diversification Rate Estimation (Höhna et al. 2019)

`dnCDBDP` can also be used to estimate the number and placement of events 
at which the rates of speciation and extinction shift from one category
to another, without the assumption that such rate shifts correspond to state
transitions in an observed character. While the number of shifts is inferred,
the number of categories N has to be specified by the user, and relates
to the length of the `lambda` and `mu` vectors. If only one of the two rates
is allowed to shift, both vectors have to be of length N (with one containing
N distinct and the other N identical elements); if both rates are allowed to
shift independently of each other, the vectors have to be of length N^2. The
`Q` matrix gives the rates of shifts between any two individual categories.
In this case, a draw from the distribution is only fixed to a tree using 
`.clamp()`.
## authors
Sebastian Hoehna
## see_also
dnCBDSP
fnCladogeneticProbabilityMatrix
fnCladogeneticSpeciationRateMatrix
## example
    # setup for a two-state ClaSSE model
    num_states = 2          # 0 and 1 are the only states

    # set basic process parameters
    root_age ~ dnUniform(0, 2)
    rho <- Probability(1/2) # sampling one half of extant lineages

    # specify extinction probabilities for each state
    mu_vec <- rep(0.1, 2)

    # Set up cladogenetic events and probabilities. Each element
    # in clado_events describes a state pattern at the cladogenetic
    # event: for example, [0, 0, 1] denotes a parent having state 0
    # and its two children having states 0 and 1.
    clado_events = [[0, 0, 1], [0, 1, 0], [1, 0, 1], [1, 1, 0]]

    # set probabilities for each cladogenetic event described above
    clado_prob <- rep(1/4, 4)

    # create cladogenetic rate matrix
    clado_matrix = fnCladogeneticProbabilityMatrix(clado_events, clado_prob,
                                                   num_states)

    # set up Q-matrix to specify rates of state changes along branches 
    q_matrix <- matrix([[0, .2],
                        [.2, 0]])

    # create vector of state frequencies at the root
    pi <- simplex([1, 2])

    # basic use of the function
    timetree ~ dnCDBDP(rootAge   = root_age,
                       lambda    = clado_matrix,
                       mu        = mu_vec,
                       Q         = q_matrix,
                       pi        = pi,
                       rho       = rho,
                       condition = "time")

## references
- citation: Barido-Sottani J, Vaughan TG, Stadler T (2020). A multitype birth–death model for Bayesian inference of lineage-specific birth and death rates. Systematic Biology, 69(5):973–986.
  doi: 10.1093/sysbio/syaa016
  url: https://academic.oup.com/sysbio/article/69/5/973/5762626
- citation: FitzJohn RG (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution, 3(6):1084-1092.
  doi: 10.1111/j.2041-210X.2012.00234.x
  url: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2012.00234.x
- citation: Goldberg EE, Igić B (2012). Tempo and mode in plant breeding system evolution. Evolution, 66(12):3701-3709.
  doi: 10.1111/j.1558-5646.2012.01730.x
  url: https://academic.oup.com/evolut/article/66/12/3701/6851227
- citation: Höhna S, Freyman WA, Nolen Z, Huelsenbeck JP, May MR, Moore BR (2019). A Bayesian approach for estimating branch-specific speciation and extinction rates. bioRxiv.
  doi: 10.1101/555805
  url: https://www.biorxiv.org/content/10.1101/555805v1.full
- citation: Maddison WP, Midford PE, Otto SP (2007). Estimating a binary character's effect on speciation and extinction. Systematic Biology, 56(5):701-710.
  doi: 10.1080/10635150701607033
  url: https://academic.oup.com/sysbio/article/56/5/701/1694265
