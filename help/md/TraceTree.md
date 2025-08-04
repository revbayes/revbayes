## name
TraceTree
## title
Trace of trees.
## description
Stores a tree trace, usually produced by passing a variable of type `Tree` to
the `mnFile` monitor in an MCMC or MCMCMC run.
## details
Important methods include:

- `getTree(n)`: returns the `n`th entry of the `TraceTree` object.
- `getTrees()`: returns a vector of trees from the `TraceTree` object, after
      excluding the burnin.
- `getUniqueTrees(alpha)`: constructs an alpha credible set of trees, i.e.,
      orders all trees in the trace by their posterior probability and adds
      them to the set one by one until the included trees accumulate fraction
      alpha of the probability mass (Huelsenbeck & Rannala 2004). By default,
      the 95% credible set is constructed (alpha = 0.95). To extract all unique
      trees from a given trace, set the argument to 1.

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
