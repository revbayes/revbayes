## name
fnHKY

## title
The Hasegawa-Kishino-Yano (1985) nucleotide rate matrix

## description
DNA evolution model proposed in Hasegawa, Kishino, and Yano (1985).

## details
In this model, nucleotides have different stationary frequencies, and transition and transversion rates are assumed to be different. Its first parameter, kappa, codes for the ratio between the rate of transitions and transversions. Its second parameter, baseFrequencies, codes for the frequencies of each nucleotide.

## authors

## see_also
fnJC
fnK80
fnK81
fnF81
fnT92
fnTrN
fnGTR

## example
	# the ratio between rates of transitions and transversions
	kappa ~ dnLognormal(0,1)
    
	# the base frequencies    
	pi ~ dnDirichlet( v(1,1,1,1) )

	# create an HKY rate matrix
	Q := fnHKY(kappa,pi)

## references
Hasegawa, M. et al. (1985). "Dating of the human-ape splitting by a molecular clock of mitochondrial DNA". Journal of molecular evolution. 22(2):160-174.
  doi: https://doi.org/10.1007/BF02101694
  url: https://link.springer.com/article/10.1007%2FBF02101694
