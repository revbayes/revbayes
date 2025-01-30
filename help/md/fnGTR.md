## name
fnGTR

## title
The General Time-Reversible rate matrix

## description
DNA evolution model proposed in Tavare (1986).

## details
In this model, states are allowed to have different stationary frequencies, and exchangeability rates between states are allowed to be different. Its first argument, exchangeRates, codes for the transition rates between states (as in other models, transition rates are assumed to be symmetric). Its second argument, baseFrequencies, codes for the stationary frequencies of these states. Note that for n states, exchangeRates should have length n*(n-1)/2, and baseFrequencies should have length n. While this is usually used for DNA (and therefore has four states), the function can take any number of states, and therefore be used for many other applications (such as aminoacid or morphological evolution).

The general time-reversible rate matrix elements will be of the form:
	Q[i, j] = c * exchangeRates[i, j] * baseFrequencies[j]

where c is a constant needed to normalize the average rate to 1.

## authors
## see_also
fnJC
fnK80
fnK81
fnT92
fnHKY
fnTrN

## example
	# exchange rates
	er ~ dnDirichlet( v(1,1,1,1,1,1) )

	# base frequencies
	pi ~ dnDirichlet( v(1,1,1,1) )
        
	# create a GTR rate matrix
	Q := fnGTR(er,pi)

## references
- citation: Tavare, S (1986). "Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences".  Lectures on Mathematics in the Life Sciences. 17:57-86
  url: http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T86.pdf
  doi: null
