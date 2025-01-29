## name
fnT92
## title
The Tamura (1992) nucleotide rate matrix
## description
DNA evolution model proposed in Tamura (1992).
## details
In this model, A and T have an equal stationary frequency, with G and C frequencies distinct, and transition and transversion rates are assumed to be different. Its first parameter, kappa, codes for the ratio between the rate of transitions and transversions. Its second parameter, gc, codes for the compound frequency of G and C nucleotides.
## authors
## see_also
fnJC
fnF81
fnK80
fnK81
fnHKY
fnTrN
fnGTR
## example
	# the ratio between rates of transitions and transversions
	kappa ~ dnExp(0.5)

	# the frequency of G and C nucleotides
	gc ~ dnUnif(0, 1)

	# create a T92 rate matrix
	Q := fnT92(kappa, gc)
## references
Tamura K (1992). "Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C-content biases". Molecular Biology and Evolution. 9:678–87.
