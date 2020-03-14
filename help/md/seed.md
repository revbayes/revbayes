## name
seed
## title
Seed set function
## description
Sets the random number generator seed given a natural number.
## details
## authors
## see_also

## example
    # pick some definitely random number
    seed(80797980)
    a <- rUniform(1,0.6,1.2)
    a
    seed(80797980)
    a <- rUniform(1,0.6,1.2)
    a # this will be the same as above!
    
## references
