## name
mvRateAgeBetaShift
## title
The RateAgeBetaShift move
## description
Resample a single node age and adjust neighboring rates to preserve distances
## details
This move first selects a tree node that is not a tip or the root of the tree.

The age of the tree node is resampled from the interval
   [max(child1.age, child2.age), parent.age]
using a Beta distribution(a,b).

The rates of the parent edge and two child edges are then modified to ensure that the rate*time
remains unchanged for the tree branches.

## authors
## see_also
mvRateAgeProposal
mvRateAgeSubtreeProposal
## example
moves.append( mvRateAgeBetaShift(tree=timetree, rates=branch_rates, tune=true, weights=n_taxa ) )
## references
