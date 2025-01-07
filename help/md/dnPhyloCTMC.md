## name
`dnPhyloCTMC`: Distribution of a phylogenetic continuous-time Markov chain

## title
The parameters of a phylogenetic model – a tree topology with branch lengths, a substitution model that describes how observations evolve over the tree, etc. – collectively form a distribution called the _phylogenetic continuous-time Markov chain_.

## description
dnPhyloCTMC gives the probability distribution of tip labels on a phylogenetic tree given an phylogenetic continuous-time Markov chain model.

## details

The likelihood of observed tip labels (specified via a clamped `AbstractHomologousDiscreteCharacterData` object) is computed using Felsenstein's pruning algorithm, with partial likelihoods stored for each branch of the tree. It is automatically outputted in the `Likelihood` column of the `mnFile()` and `mnScreen()` monitors (which can be suppressed with `likelihood = FALSE`).

## authors
## see_also
- Tutorial on [graphical models](https://revbayes.github.io/tutorials/intro/graph_models)

- Tutorial on [specifying a phylogenetic continuous-time Markov chain](https://revbayes.github.io/tutorials/ctmc/) model

## example

```rb
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
```

## references
- citation: Felsenstein J., 1973. Maximum Likelihood and Minimum-Steps Methods for Estimating Evolutionary Trees from Data on Discrete Characters. Systematic Biology 22:3, 240--249
  doi = 10.1093/sysbio/22.3.240
- citation: Felsenstein, J. (1981). Evolutionary trees from DNA sequences: A maximum likelihood approach. Journal of Molecular Evolution. 17 (6): 368–376.
  doi : 10.1007/BF01734359
- citation: Höhna, S., Landis, M.J. and Heath, T.A. 2017. Phylogenetic inference using `RevBayes`. Curr. Protoc. Bioinform.
57:6.16.1-6.16.34.
  doi: 10.1002/cpbi.22
  url: null


