## name
dnFossilRecord
## title
## description
The fossil record of a set of species, conditioned on a fossilized birth-death range process: the probability of the observed fossil occurrences given the species ranges that produced them.
## details
This is the observation half of the fossilized birth-death range model. The range process (`dnFBDRP` or `dnFBDSP`) describes how species ranges diversify; this distribution gives the probability of the fossil record itself.

The `complete` argument selects how sampled specimens make it into the reported record: `complete=TRUE` (every sampled occurrence is reported) or `complete=FALSE` (only the oldest and youngest occurrence of each species are reported, i.e. first/last). It also determines the support of the augmented oldest age, which belongs to the range process, so both nodes use the model given here.

The occurrences, the fossil sampling rate, its timeline, and the augmented occurrence ages are all read from the range process rather than given as arguments here. The sampling rate also appears in the range process's non-detection term, so the two nodes must share it.

The deprecated `dnFBDRMatrix` fuses this term and the range process into a single node, and additionally offers the truncated (exchangeable-occurrence) reporting model via its `truncated=K` cap.
## authors
June Walker
## see_also
dnFossilizedBirthDeathRange
dnBDS
mvResampleFBDR
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

# the birth-death range process
bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, clamped to the observed occurrences
rec ~ dnFossilRecord(ranges=bd, complete=false)
rec.clamp(taxa)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
