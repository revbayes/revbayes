## name
readTreeTrace
## title
Function to read in a tree trace, usually produced as the output of an MCMC.
## description
Reads trees (Nexus or Newick accepted) from a file or folder containing a set of trees and saves them in one object. 
## details
Either a file name or a directory must be provided as argument. If a folder is provided, all the files that contain trees in that directory are read in the same object.

`tree_trace = readTreeTrace(..., nruns = 1)` returns a `TreeTrace` object.
`tree_trace = readTreeTrace(..., nruns > 1)` returns a `TreeTrace[]` object, in which `trace[n]` corresponds to the trace of run $n$.

A `TreeTrace` stores every `thinning`th sample from a file, starting at the sample numbered `offset + 1`.

Pertinent methods of a `TreeTrace` object include:
`tree_trace.getTree(n)`: returns the $n$th entry of the `TreeTrace` object.
`tree_trace.getTrees()`: returns a vector of trees from the `TreeTrace` object, after excluding the burnin.


## authors
## see_also
readTrace
readCharacterData
readTrees
## example
    # read a tree trace
    tree_trace = readTreeTrace("my_filename.tree", treetype="clock", burnin=0.5)

    # Ignore the first 10 samples
    thinned_trees = readTreeTrace("my_filename.tree", offset = 10, thinning = 10, burnin = 0.5)
    
    thinned_trees.getTree(1) # Returns the 11th tree (offset + 1) in the file
    thinned_trees.getTree(2) # Returns the 21st tree (offset + 1 + thinning) in the file

    thinned_trees.getTrees()[1] # Returns the first sampled tree after excluding the burnin fraction

    # make a summary MCC tree
    mcc_tree = mccTree(trace=tree_trace, file="mcc.tree")

## references
