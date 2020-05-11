## name
dnMultivariateNormal
## title
Multivariate Normal Distribution
## description
The multivariate normal distribution has the probability density:

f(x) = det(2 pi Sigma)^(-1/2) e^{-(1/2) (x-mu)' Sigma^-1 (x-mu)}

where mu is a vector of mean values and Sigma is a covariance matrix. Note, this distribution may also be parameterized in terms of the precision matrix, Sigma^-1.

## details
## authors
Sebastian Hoehna
Michael Landis
## see_also
dnNormal
dnWishart
mvCorrelationMatrixUpdate
## example
	dim = 4
	df = 100
	kappa <- 2
	Sigma ~ dnWishart(df, kappa, dim)
	for (i in 1:dim) { mu[i] ~ dnUnif(-1, 1) }
	x ~ dnMultivariateNormal( mean=mu, covariance=Sigma )
	mv[1] = mvCorrelationMatrixElementSwap(Sigma)
	mv[2] = mvCorrelationMatrixRandomWalk(Sigma)
	mv[3] = mvCorrelationMatrixSingleElementBeta(Sigma)
	mv[4] = mvCorrelationMatrixSpecificElementBeta(Sigma)
	mv[5] = mvCorrelationMatrixUpdate(Sigma)
	mv[6] = mvVectorSlide(x)
	
## references
