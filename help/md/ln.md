## name
ln
## title
Natural log function
## description
Returns the natural log of a (positive) value.
## details
## authors
Jordan Koch
## see_also
log
## example
	# create a stochastic node with an exponential distribution
	x ~ dnExponential(1)
	
	# create a determinstic node that takes the natural log of x
	y := ln(x)
	
	# print the values for x and y
	x # x has the stochastic value of 2.940149
	y # y has the determined value of 1.07846
	
## references
