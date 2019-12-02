## name
v
## title
Create a vector
## description
'v' creates a vector of the elements '...'
## details
'v' creates a vector of the elements '...', which are objects of a common base type. Vector elements may themselves be vectors.
## authors
Michael Landis
## see_also
simplex
rep
## example
	# create a vector, Natural[]
	x <- v(1,2,3)
	x <- x + 1
	x
	
	y <- v(2,4,6)
	# create a vector of Natural[][]
	z <- v(x,y)
	z
	z[0]
	
## references
