## name
dnPhyloCTMC
## title
Distribution of a phylogenetic continuous-time Markov chain
## description
Gives the probability distribution of the character state vectors at the leaves
of a phylogenetic tree, given a phylogenetic continuous-time Markov chain
model.
## details
The parameters of a phylogenetic model -- a tree topology with branch lengths,
a substitution model that describes how observations evolve over the tree, etc.
-- collectively form a distribution called the _phylogenetic continuous-time
Markov chain_.

The likelihood of observed character state vectors (specified via clamping the
distribution to a `AbstractHomologousDiscreteCharacterData` object) is computed 
using Felsenstein's pruning algorithm, with partial likelihoods stored for each
branch of the tree. It is automatically outputted in the `Likelihood` column of
the `mnFile()` and `mnScreen()` monitors (which can be suppressed with
`likelihood = FALSE`).

For more details, see the tutorials on [graphical models](https://revbayes.github.io/tutorials/intro/graph_models) and on 
[specifying a phylogenetic continuous-time Markov chain](https://revbayes.github.io/tutorials/ctmc/) model.
## authors
## see_also
## example
    # Read character data from a file
    chars <- readDiscreteCharacterData("myData.nex")
    taxa = chars.taxa()
    
    # Draw a tree with branch lengths
    tree ~ dnUniformTopologyBranchLength( taxa, branchLengthDistribution=dnExp(10.0) )
    
    # Define a rate matrix
    q_matrix <- fnJC(4)
    
    # Create stochastic node with the tip distribution given by `tree` and `q_matrix`
    x ~ dnPhyloCTMC(tree = tree, Q = q_matrix)
    
    # Clamp observed characters to the node
    x.clamp(chars)
    
    # Calculate the probability of the observed characters under the given distribution
    x.lnProbability()

    # Simulate characters
    sim ~ dnPhyloCTMC(tree = tree, Q = q_matrix, nSites = 24)
    
    # Print simulated characters to screen
    sim.show()
    
    # Write dataset to file
    writeNexus("simulatedData.nex", sim)

## references
- citation: Felsenstein J (1973). Maximum likelihood and minimum-steps methods for estimating evolutionary trees from data on discrete characters. Systematic Biology, 22(3):240--249.
  doi: 10.1093/sysbio/22.3.240
- citation: Felsenstein J (1981). Evolutionary trees from DNA sequences: A maximum likelihood approach. Journal of Molecular Evolution, 17(6):368--376.
  doi: 10.1007/BF01734359
- citation: Hoehna S, Landis MJ, Heath TA (2017). Phylogenetic inference using `RevBayes`. Current Protocols in Bioinformatics, 57:6.16.1--6.16.34.
  doi: 10.1002/cpbi.22
