## name
fnFreeK
## title
Free K Rate Matrix
## description
This function generates and returns a free rates matrix. 
## details
This function accepts both a vector or a matrix of non-negative, real numbers in the 
first argument to automatically generate a rate matrix with corresponding substitution 
rates, returning a rate matrix object. The function will fill rates in the matrix from
left to right as provided in the first argument, skipping the diagonal when using a vector
as input. For this reason,using a vector of lengths 2 to 5 will create a 2-by-2 rate
matrix but a vector of length 6 will create a 3-by-3 rate matrix as fnFreeK will have 
enough values to fill the matrix. Using a matrix to create in fnFreeK will create a 
rate matrix object with rates filled in their respective position in the provided matrix.
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
    # Define vector to pass, this will create a 2-by-2 matrix
    x <- v(0.5, 0,5)
    # Use fnFreeK to create rate matrix
    # Note the second argument is true in this case so rates will be normalized
    fnFreeK(x)
    [ [ -1.0000, 1.0000 ] ,
      [ 1.0000, -1.0000 ] ]
    
    # Case where rates are not normalized
    x <- v(0.5, 0.5)
    fnFreeK(x, false)
    [ [ -0.5000, 0.5000 ] ,
      [ 0.5000, -0.5000 ] ]

    # Define matrix for 3-by-3 rate matrix object
    x <- v([0, .6, .4], [.2, 0, .4], [.3, .3, 0])
    # Create rate matrix object
    fnFreeK(x)
    [ [ -1.0000, 0.6000, 0.4000 ] ,
      [ 0.2000, -0.6000, 0.4000 ] ,
      [ 0.3000, 0.3000, -0.6000 ] ]
## references
