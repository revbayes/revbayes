## name
structure
## title
The structure of a variable
## description
Shows all the information about a given variable.
## details
## authors
Sebastian Hoehna
## see_also
type
## example
	# create a variable
	a <- 1
	b := exp(a)
	# now create a deterministic variable which will be a child of 'b'
	c := ln(b)
	# now create a constant variable which will not be a child of 'b'
	d <- ln(b)
	
	str(b)
	
## references
