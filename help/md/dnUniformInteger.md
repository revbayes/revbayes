## name
dnUniformInteger
## title
Uniform Integer Distribution
## description
This distribution creates a stochastic node drawing a random integer from a uniform distribution.
## details
This distribution will randomly draw an integer using uniform distribution
from a minimum to maximum integer set in the first and second arguments.
This function can alsonbe called using the alias 'dnUnifInt'.
## authors
Sebastion Hoehna
## see_also
dnNormal
dnExponential
## example
    # Create and assign stochastic node
    # To obtain a new value for x the distribution will need to be called again
    x ~ dnUniformInteger(1, 10)
## references
