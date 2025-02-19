## name
dnUniformInteger
## title
Uniform Integer Distribution
## description
This function creates a stochastic node drawing a random integer from a uniform distribution.
## details
This function will randomly draw an integer from a uniform distribution
from a minimum and maximum integer set in the first and second arguments.
dnUniformInteger must be defined as a stochastic node. This function can also
be called using the alias 'dnUnifInt'.
## authors
Sebastion Hoehna
## see_also
dnNormal
dnExponential
## example
    # Create stochastic node
    x ~ dnUniformInteger(1, 10)
    # See what x was assigned
    x
    5
## references
