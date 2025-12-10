## name
dnPhyloCTMCDASiteIID
## title
Data augmentation-based phylogenetic continuous-time Markov chain distribution with IID sites
## description
The probability distribution of the character state vectors at the tips of a phylogenetic tree, given a phylogenetic continuous-time Markov chain model.
## details
The parameters of a phylogenetic model -- a tree topology with branch lengths, a substitution model that describes how observations evolve over the tree, etc. -- collectively form a distribution called the _phylogenetic continuous-time
Markov chain_. In this distribution, the branch history is sampled (through data augmentation), instead of integrated out (as in `dnPhyloCTMC`). An instance of character history can therefore be retrieved from the distribution (using `.characterHistories`) The probability of observed character state vectors (specified via clamping the
distribution to a `AbstractHomologousDiscreteCharacterData` object) is computed by summing the probability of the history on each branch. The probability of the history on a branch is the product of the probabilities of waiting times between events (or the probability of no event in the final segment) given the current rate of change (see May and Moore 2020).

Note that when `rootFrequencies` is not provided, the distribution assumes stationary frequencies at the root.
The stationary frequencies will be calculated numerically if the tree has a root branch, and analytically otherwise.
When `rootFrequencies` is provided, then stationarity at the root will not be assumed, and the likelihood calculation will be based on the provided root frequencies.

If this distribution is used together with `dnPhyloOUSD` or `dnPhyloBMSD`, specify `nSites=1`.

## authors
Priscilla Lau
## see_also
mvCharacterHistory
## example
    # Read character data from a file
    chars <- readDiscreteCharacterData("myData.nex")
    taxa = chars.taxa()

    # Draw a tree with branch lengths
    tree ~ dnUniformTopologyBranchLength( taxa, branchLengthDistribution=dnExp(10.0) )

    # Define a rate matrix
    q_matrix <- fnJC(4)

    rf <- rep(1/4, 4)
    
    # Create stochastic node with the tip distribution given by `tree` and `q_matrix`
    x ~ dnPhyloCTMCDASiteIID(tree = tree, Q = q_matrix, rootfrequencies=rf)
    char_hist := x.characterHistories()

    # Clamp observed characters to the node
    x.clamp(chars)

    # Calculate the probability of the observed characters under the given distribution
    x.lnProbability()

## references
- citation: May MR, Moore BR (2020). A Bayesian approach for inferring the impact of a discrete character on rates of continuous-character evolution in the presence of background-rate variation. Systematic biology, 69(3):530-544.
  doi: 10.1093/sysbio/syz069
  url: https://academic.oup.com/sysbio/article/69/3/530/5609130
