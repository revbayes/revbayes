## name
dnLKJ
## title
LKJ Distribution
## description
The LKJ (Lewandowski-Kurowicka-Joe) distribution on correlation matrices with concentration parameter eta.
## details
The LKJ distribution is uniform over positive-definite correlation matrices when eta=1.The probability density of a correlation matrix under the LKJ distribution is:f(x) = det(x)^(eta - 1)
## authors
Michael R. May
## see_also
dnLKJPartial
## example
	
	# we simulate a correlation matrix.
	R <- rLKJ(n=1, eta=1, dim=5)
	
	# let's print the simulated correlation matrix
	R
	
## references
- citation: Lewandowski D, D Kurowicka, H Joe (2009). Generating random correlation
    matrices based on vines and extended onion method. Journal of Multivariate Analysis,
    100(9):1989-2001.
  doi: null
  url: null
