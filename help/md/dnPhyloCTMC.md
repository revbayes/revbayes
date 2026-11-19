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

Optionally, an observation error model can be applied to 
account for scoring ambiguity (e.g., in morphological datasets). 
This distinguishes between the true biological state and the recorded score.
When `observationErrorProbability` (epsilon) is > 0, the tip likelihoods
are initialized as a mixture: with probability (1 - epsilon), the score 
is accurate; with probability epsilon, the score is drawn from the 
distribution defined by `observationErrorFrequencies`

The `coding` argument specifies which site patterns could have been sampled and corrects the likelihood for excluded patterns (ascertainment bias correction). The default, `"all"`, applies no correction. For Standard data, you can use `"all"`, `"informative"`, `"variable"`, or `"informative|ordered"`, but only one at a time. Other data types accept only `"all"`, except for Binary and Restriction data.

For `type="Binary"` or `type="Restriction"`, you can combine options in a string separated by `"|"` or as a vector. For example, `coding="noabsencesites|nosingletonpresence"` and `coding=["noabsencesites", "nosingletonpresence"]` both exclude sites where all taxa have 0 or exactly one taxon has 1. Each option may appear only once. The options `"variable"` and `"nosingletons"` each combine the corresponding absence and presence exclusions, and `"informative"` combines `"variable"` with `"nosingletons"`. Combining codings for other data types gives an error.

For Standard data, `coding="informative"` computes the exact probability of all parsimony-uninformative site patterns for any number of states using a dynamic programming algorithm related to that of Koch and Holder (2012). Memory use grows exponentially with the number of states, more than doubling with each additional state. Characters with many states, for example 10, can require many gigabytes of memory and cause the analysis to fail at the start. Memory use grows linearly with the number of taxa, the number of site rate categories or other site mixture components, and the number of distinct patterns of missing data across characters.

For Standard data, you can add `"ordered"` to `"informative"` using `coding="informative|ordered"` or `coding=["informative", "ordered"]`. Plain `"informative"` excludes patterns that are uninformative under unordered parsimony. With `"ordered"`, the correction instead excludes patterns that are uninformative under ordered (Wagner) parsimony. Under ordered parsimony, a site is informative if there is some threshold between two adjacent states with at least two taxa on each side. A taxon with an ambiguous score that spans the threshold counts on neither side. Use this option when your characters are ordered, for example with `Q=fnOrderedRateMatrix(5)`. Using `"ordered"` on its own or with any coding other than `"informative"` gives an error.

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
- citation: Höhna S, Landis MJ, Heath TA (2017). Phylogenetic inference using `RevBayes`. Current Protocols in Bioinformatics, 57:6.16.1--6.16.34.
  doi: 10.1002/cpbi.22
- citation: Koch JM, Holder MT (2012). An algorithm for calculating the probability of classes of data patterns on a genealogy. PLoS Currents, 4:e4fd1286980c08.
  doi: 10.1371/4fd1286980c08
