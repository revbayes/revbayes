## name
mvVectorSlide
## title
Additive Move on a Vector
## description
The mvVectorSlide move applies a sliding-window proposal to all elements of a vector of continuous parameters simultaneously.
It works by drawing a single random number from a uniform distribution and adding that same number to every element of the vector
## details
This move has two main arguments:delta — controls how large the proposed changes can be (the size of the sliding window).
weight — determines how often this move is selected during MCMC. This means all elements shift up or down by the same amount, keeping their relative differences constan. Useful when you want to adjust the overall level of a vector while preserving the relative differences among its elements
## authors
## see_also
mvVectorScale
mvSlide
## example
## references

