## name
fnFreeK
## title
Free K Rate Matrix
## description
This function generates and returns a free rates matrix. 
## details
This function accepts both RealPos[] or RealPos[][] as the first argument to automatically
generate a rate matrix with corresponding substitution rates, returning a rate matrix object.
Users can specify if matrix should be normalized in the second argument using a boolean 
variable (default TRUE). Lastly users can specify what matrix exponential method to
use (default eigen) with a string. Possible options include:
scalingAndSquaring
scalingAndSquaringPade
scalingAndSquaringTaylor
uniformization
eigen
## authors
Michael Landis
## see_also
RateMatrix
fnFreeBinary
fnFreeSymmetricRateMatrix
## example
    # Define vector to pass
    x <- [1, 1, 1, 1]
    # Use fnFreeK to create rate matrix
    fnFreeK(x)
    [ [ -1.0000, 1.0000 ] ,
      [ 1.0000, -1.0000 ] ]
## references
