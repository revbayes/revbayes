## name
quantile
## title
Quantile function
## description
Calculates the sample quantiles corresponding to the given probability.
The smallest observation corresponds to a probability of 0.0 and the 
largest to a probability of 1.0. The median corresponds to a probability of
0.5.
## details
## authors
## see_also
`median`
## example
    b <- v(2,4,6,7,9)
    quantile(b, k = 0.0)
    # returns 2
    quantile(b, k = 0.5)
    # returns 6
    quantile(b, k = 1.0)
    # returns 9
    
## references
