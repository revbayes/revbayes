## name
TraceTree
## title
Trace of trees.
## description
Stores a tree trace, usually produced by passing a variable of type `Tree` to
the `mnFile` monitor in an MCMC or MCMCMC run.
## details
Individual trees can be accessed as follows:

- `getTree(n)`: returns the `n`th tree.
- `getTree(n, post=TRUE)`: returns the `n`th post-burnin tree.
- `getTrees()`: returns the entire vector of post-burnin trees.
      
The methods `.computeEntropy()`, `.computePairwiseRFDistances`,
`.getUniqueTrees()`, `.isTreeCovered()`, and `.summarize()` first construct
a credible set of trees, and then perform the corresponding operation on this
set. A credible set is constructed by ordering all trees in the trace by their
posterior probability, and adding them to the set one by one until the included
trees accumulate the desired probability specified by `credibleTreeSetSize`.
By default, this argument is set to 0.95; to perform a given operation on all
trees in the trace, set it to 1.

The credible set methods also take a `probabilistic` argument that determines
how the credible set should be constructed when it is not possible to obtain a
cumulative probability precisely equal to `credibleTreeSetSize`. This will be
the case when the probability p of the first n trees falls short of it, and the
probability (p + q) of the first n + 1 trees exceeds it (Huelsenbeck & Rannala
2004). If `probabilistic=TRUE` (default), the (n + 1)-th tree will be included
in the set with probability (`credibleTreeSetSize` - p) / q. It will always be
included if `probabilistic=FALSE`. Note that when `probabilistic=TRUE`, the
methods in question may yield different results under different random seeds,
and they may even construct an empty credible set a certain fraction of the
time if the trace is dominated by a particular tree whose probability exceeds
`credibleTreeSetSize`. An exception is thrown when an empty credible set is
encountered by a method that does not expect one.
## authors
## see_also
mnFile
readAncestralStateTreeTrace
readTreeTrace
Trace
Tree
## example
    # read a tree trace and ignore the first 10 samples
    thinned_trees = readTreeTrace("my_filename.tree", offset = 10, thinning = 10, burnin = 0.5)
    
    thinned_trees.getTree(1) # Returns the 11th tree (offset + 1) in the file
    thinned_trees.getTree(2) # Returns the 21st tree (offset + 1 + thinning) in the file

    thinned_trees.getTrees()[1] # Returns the first sampled tree after excluding the burnin fraction

## references
- citation: Huelsenbeck JP, Rannala B (2004). Frequentist properties of Bayesian posterior probabilities of phylogenetic trees under simple and complex substitution models. Systematic Biology, 53(6):904-913.
  doi: 10.1080/10635150490522629
  url: https://academic.oup.com/sysbio/article-abstract/53/6/904/1651356
