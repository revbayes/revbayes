## name
fnFreeBinary

## title
Free Binary transition rate matrix

## description
Constructs a transition rate matrix between two states.

## details
This function enables the user to define the non-normalized off-diagonal elements of the matrix while also providing the option to normalize the matrix. Normalization ensures that the mean instantaneous rate equals 1. For example, a branch of length 1 in a non-clock tree corresponds to an expected 1 substitution per character.
 
 It takes in two arguments: 
    (1) transition_rates (tr) - A vector of real numbers of length 2 that represents the rate of transition between states.
    (2) rescaled - A boolean value that indicates whether or not the matrix should be normalized. Takes on TRUE by default.

0 and 1 represent our 2 states:
Q = [[-q_{01}, q_{01}],[q_{10}, -q_{10}]]

## authors

## see_also
fnFreeK

## example
# Under the ERM model
rate_pr := phylogeny.treeLength() / 10
mu ~ dnExp(rate_pr)

moves.append( mvScale(mu, lambda=1, weight=2.0) )

rate := rep(mu, 2)

Q := fnFreeBinary(rate, rescaled=false)

## references