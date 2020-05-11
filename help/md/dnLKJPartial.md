## name
dnLKJPartial
## title
LKJ Distribution (for partial correlation matrices)
## description
The LKJ (Lewandowski-Kurowicka-Joe) distribution (on the partial correlation matrix) with concentration parameter eta.
## details
The LKJPartial distribution is uniform over positive-definite correlation matrices when eta=1.

The probability density of a correlation matrix under the LKJ distribution is:

f(x) = det(x)^(eta - 1)

## authors
Michael R. May
## see_also
dnLKJ
## example
	# we simulate a partial correlation matrix.
	P <- rLKJPartial(n=1, eta=1, dim=5)
	
	# let's print the simulated partial correlation matrix
	P
	
## references
- citation: Lewandowski D, D Kurowicka, H Joe (2009). Generating random correlation
    matrices based on vines and extended onion method. Journal of Multivariate Analysis,
    100(9):1989-2001.
  doi: null
  url: null
