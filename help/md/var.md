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
\[
s^2 = \frac{1}{n-1} \sum_{i=1}^{n} (x_i - \bar{x})^2
\]
Where:
- \( s^2 \) = sample variance  
- \( n \) = number of data points  
- \( x_i \) = each individual data point  
- \( \bar{x} \) = mean of the data 

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
