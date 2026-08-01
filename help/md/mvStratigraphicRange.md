## name
mvStratigraphicRange
## title
Stratigraphic range move
## description
Resamples the first and last appearances of one randomly chosen taxon in a fossilized birth-death speciation process (`dnFBDSP`).
## details
A stratigraphic range is bounded by its first and last appearances. Both are latent: the occurrence record reports the bins that contain them, not the ages themselves. A tree carries divergence times only, so these two ages are state of the distribution rather than elements of its value, and no move on the tree can reach them.

This move redraws both, independently and uniformly within the bins the record reports. The support is fixed by the data rather than by the current state, so the proposal is an independence sampler and carries no Hastings ratio.

The move applies to `dnFBDSP` only. A `dnFBDRP` keeps the same two ages in columns 2 and 3 of its value, where the generic matrix element moves sample them.

Without this move the appearances stay at their initial draw and the chain samples the wrong distribution, so the process warns at startup when no such move is attached.
## authors
June Walker
## see_also
dnFossilizedBirthDeathSpeciation
dnFossilRecord
mvRotateNode
## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, taxa=taxa)

rec ~ dnFossilRecord(ranges=tr, complete=false)
rec.clamp(taxa)

moves.append( mvStratigraphicRange(tr, weight=taxa.size()) )
## references
	- citation: The fossilized birth-death model for the analysis of stratigraphic range data under different speciation modes. Stadler, Tanja et al. Journal of theoretical biology, 447:41-55.
	  doi: doi.org/10.1016/j.jtbi.2018.03.005
	  url: https://www.sciencedirect.com/science/article/pii/S002251931830119X
