## name
sqrt
## title
Square root of a number
## description
The 'sqrt' function return the square root of a number.
## description
Takes the square root of some positive number `x`.
## details
## authors
## see_also
`power`
## example
        # compute the square root of a real number
        x <- 3.0
        root <- sqrt(x)
        if ( abs(root*root - x) > 1.0e-15) {
            print("Problem computing the square root.")
        } else {
            print("Correct computation of the square root.")
        }
## references
