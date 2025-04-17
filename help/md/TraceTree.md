## name
TraceTree
## title
Trace of trees.
## description
Stores a tree trace, usually produced by passing a variable of type `Tree` to the `mnFile` monitor in an MCMC or MCMCMC run.
## details
Important methods include:

- `getTree(n)`: returns the `n`th entry of the `TraceTree` object.
- `getTrees()`: returns a vector of trees from the `TraceTree` object, after excluding the burnin.

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
