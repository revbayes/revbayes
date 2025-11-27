## name
var
## title
Variance
## description
Calculate variance of a vector of real numbers
## details
This function accepts a vector of real numbers and returns the variance.
This a measure of how the data points deviate from the mean which is calculated
as follows:

s^2 = (1 / (n - 1)) * \sum (x[i] - \bar{x})^2

which is essentially:
(sum of squared differences from the mean) / (n - 1)

## authors
## see_also
mean
stdev
median
## example
    # Define vector to calculate variance
    x <- v(1, 2, 3, 4)
    # or
    x <- [1, 2, 3, 4]
    # Calculate variance
    var(x)
## references
