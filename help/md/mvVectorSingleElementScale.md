## name
mvVectorSingleElementScale
## title
Scale Move for a Single Element of a Vector
## description
mvVectorSingleElementScale proposes a random scaling of one element from a vector of positive values
## details
This move picks one element from a vector and multiplies it by a random factor. 
The factor is drawn from:scale = exp(lambda * (u - 0.5)), where u ~ Uniform(0,1).
This move is best for vectors with positive values, but it can technically be used with negative values as well
## authors
## see_also
## example
## references
