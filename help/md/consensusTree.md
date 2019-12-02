## name
consensusTree
## title
## description
Calculates the majority-rule consensus topology from a trace of trees and summarizes branch lengths.
## details
## authors
Seraina Klopfstein
Will Freyman
Will Pett
Sebastian Hoehna
## see_also
mapTree
mccTree
treeTrace
readTreeTrace
## example
	# Read in tree trace
	tree_trace = readTreeTrace("output/my.trees", burnin=0.25)
	
	# Generate the majority-rule consensus tree
	map_tree = consensusTree(trace=tree_trace, cutoff=0.5, file="consensus.tree")
	
## references
