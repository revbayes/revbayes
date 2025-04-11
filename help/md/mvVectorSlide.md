## name
mvVectorSlide
## title
Additive Move on a Vector
## description
mvVectorScale proposes multiplicative changes to all elements of a rate vector at once
## details
This proposal randomly slides all elements of a vector using the same sliding factor
A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5) and slides the current vale by a sliding offset delta  = ( lambda * u )
where lambda is the tuning parameter of the proposal to influence the size of the proposals
## authors
## see_also
mvSlide
## example
## references
