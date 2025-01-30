## name
fnK81

# title
The Kimura (1981) nucleotide rate matrix

## description
DNA evolution model proposed in Kimura (1981).

## details
In this model, transition and transversion rates are allowed to be different, and transversion rates for A <-> C, G <-> T and A <-> T, C <-> G transversions are different as well. The first argument, kappa1, defines the ratio between the rate of transitions and the rate of A <-> C, G <-> T transversions. The second argument, kappa2, defines the ratio between the rate of A <-> T, C <-> G transversions and the rate of A <-> C, G <-> T transversions. The third argument, baseFrequencies, defines the stationary frequencies of nucleotide bases. Note that the original Kimura (1981) model assumed equal base frequencies, so this function is more general (if ran without a baseFrequencies argument, however, this is equivalent to K81, since the default is all frequencies equal). 

The K81 rate matrix elements will be of the form:
Q[i, j] = c, if i<->j is an A<->C/G<->T transversion
		= c * kappa1, if i<->j is a transition
		= c * kappa2, if i<->j is an A<->T/C<->G transversion

where c is a constant needed to normalize the average rate to 1. If using the baseFrequencies parameter, those elements are multiplied by baseFrequencies[j].

## authors
## see_also
fnJC
fnK80
fnF80
fnT92
fnHKY
fnTrN
fnGTR

## example
	# the ratio between rates of transitions and A<->C/G<->T transversions
	kappa1 ~ dnExp(0.5)

	# the ratio between rates of A<->T/C<->G and A<->C/G<->T transversions
	kappa2 ~ dnExp(0.5)

	# create a K81 rate matrix
	Q := fnK81(kappa1, kappa2)

	# base frequencies 
	baseFrequencies ~ dnDirichlet(v(1,1,1,1))

	# K81 rate matrix with non-equal base frequencies
	Q := fnK81(kappa1, kappa2, baseFrequencies)

## references
- citation: Kimura M (1981). "Estimation of evolutionary distances between homologous nucleotide sequences". Proceedings of the National Academy of Sciences of the United States of America. 78:454–8.
  doi: https://doi.org/10.1073/pnas.78.1.454
  url: https://www.pnas.org/doi/abs/10.1073/pnas.78.1.454
