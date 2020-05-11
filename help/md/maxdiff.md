## name
maxdiff
## title
## description
Finds the maximum difference in clade probabilities between two posterior samples.
## details
## authors
Will Pett
## see_also
readTreeTrace
## example
	# Read in tree trace
	tree_trace = readTreeTrace("output/my.trees", burnin=0.25, nruns=2)
	
	# Compute the maxdiff statistic
	maxdiff = maxdiff(traces=tree_trace)
	
## references
