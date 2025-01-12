## name
fnMinBLTimeScaling
## title
## description
Time-scales an undated tree based on a vector of tip ages using the minimum
branch length ("MBL") approach (Laurin 2004; Bapst 2014).
## details
The age of each internal node is based on the age of the oldest tip descended
from it. However, if t0 is the oldest tip descended from a given node (denoted
x) and also the oldest tip descended from that node's parent (denoted y), then
setting t(x) = t(y) = t(t0) would produce zero-length branches y->x and x->t0.
We avoid this by requiring every branch to be no shorter than some user-supplied
constant. This has the effect of shifting node ages deeper into the past.

Conceptually, the undated tree would usually correspond either to a bare
topology (a tree without branch lengths) or a tree with branch lengths in units
of expected change; in practice, both `BrachLengthTree` and `TimeTree` arguments
are accepted. In this implementation of the MBL approach, both terminal and
internal branches are required to be greater than or equal to the specified
minimum. If there is uncertainty associated with the age of a given tip,
the midpoint of the uncertainty range is used for time-scaling.

The algorithm is not stochastic (i.e., it always returns the same time-scaled
tree for a given input), and is primarily intended to generate a plausible
starting tree for MCMC analyses.
## authors
David Černý
Laura Mulvey
## see_also
simStartingTree
## example
    # Read in an undated tree
    undated_tree <- readTrees("undated.nex")[1]
    
    # Read tip age data from a file
    taxa <- readTaxonData("tipages.tsv")
    
    # Time-scale using a minimum branch length of 3 Myr
    dated_tree <- fnMinBLTimeScaling(undated_tree, taxa, 3.0)
    
    print(undated_tree) # The original tree remains unchanged
    print(dated_tree)   # A new, dated tree has been returned
    
## references
- citation: Bapst DW (2014). Assessing the effect of time-scaling methods on
phylogeny-based analyses in the fossil record. Paleobiology, 40(3):331-351.
  doi: 10.1666/13033
  url: null
- citation: Laurin M (2004). The evolution of body size, Cope's rule and the
origin of amniotes. Systematic Biology, 53(4):594-622.
  doi: 10.1080/10635150490445706
  url: null
