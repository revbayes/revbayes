## name
readTreeTrace
## title
Function to read in a tree trace, usually produced as the output of an MCMC.
## description
Reads trees (Nexus or Newick accepted) from a file or folder containing a set of trees and saves them in one object. 
## details
Either a file name or a directory must be provided as argument. If a folder is provided, all the files that contain trees in that directory are read in the same object.

`tree_trace = readTreeTrace(..., nruns = 1)` returns a `TraceTree` object.
`tree_trace = readTreeTrace(..., nruns > 1)` returns a `TraceTree[]` object, in which `trace[n]` corresponds to the trace of run `n`.

The resulting `TraceTree` object stores every `thinning`th sample from a file, starting at the sample numbered `offset + 1`.

## authors
## see_also
readTrace
readCharacterData
readTrees
## example
    # read a tree trace
    tree_trace = readTreeTrace("my_filename.tree", treetype="clock", burnin=0.5)

    # make a summary maximum clade credibility (MCC) tree
    mcc_tree = mccTree(trace=tree_trace, file="mcc.tree")

## references
