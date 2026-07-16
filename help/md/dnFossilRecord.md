## name
dnFossilRecord
## title
## description
The fossil record of a set of species, conditioned on a fossilized birth-death range skeleton: the probability of the observed fossil occurrences given the species ranges that produced them.
## details
This is the observation half of the fossilized birth-death range model. The skeleton (`dnFBDRP`) describes how species ranges diversify; this distribution supplies the probability of the fossil record itself, so the occurrences enter the model as clamped data rather than as a constructor argument.

The `reporting` argument selects the retention model, i.e. how sampled specimens make it into the reported record: `complete` (every sampled occurrence is reported), `firstlast` (only the oldest and youngest occurrence of each species are reported) or `uniform` (the reported occurrences are an exchangeable subset, capped at the largest observed count). The reporting model is pushed onto the skeleton, since it also determines the support of the augmented oldest age.

The fossil sampling rate psi, its timeline, and the augmented occurrence ages are read from the skeleton rather than taken as arguments here: psi also appears in the skeleton's non-detection term, so the two nodes must share it.

The deprecated `dnFBDRMatrix` fuses this term and the skeleton into a single node.
## authors
Walker Pett
## see_also
dnFossilizedBirthDeathRange
mvResampleFBDR
## example
lambda ~ dnExp(10)
mu ~ dnExp(10)
psi ~ dnExp(10)

# the birth-death range skeleton
bd ~ dnFBDRP(lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

# the fossil record, clamped to the observed occurrences
rec ~ dnFossilRecord(skeleton=bd, reporting="firstlast", taxa=taxa)
rec.clamp(taxa)

moves.append( mvMatrixElementScale(bd, weight=taxa.size()) )
moves.append( mvMatrixElementSlide(bd, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
