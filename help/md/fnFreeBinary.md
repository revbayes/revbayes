## name
fnFreeBinary

## title
Free Binary transition rate matrix

## description
Constructs a transition rate matrix between two states.

## details
 This function creates a transition rate matrix (between two states) for implementation in modeling evolutionary processes. This function accepts non-normalized transition rates which can then be normalized. This normalization process is a key component of ensuring that the total rates (rows) sum to 0 (requirement for transition rate matrices).
 
 It takes in two arguments: 
    (1) transtion_rates (tr) - A real number that represents the rate of transition between states.
    (2) rescaled - A boolean vlue that indicates whether or not the matrix should be normalized. Takes on TRUE by default.

0 and 1 represent our 2 states:
Q = [[-q_{01}, q_{01}],[q_{10}, -q_{10}]]

## authors
Sigournie Brock

## see_also
fnFreeK

## example
# Under the ERM model
rate_pr := phylogeny.treeLength() / 10
mu ~ dnExp(rate_pr)

moves.append( mvScale(mu, lambda=1, weight=2.0) )

NUM_STATES = 2
NUM_RATES = NUM_STATES * (NUM_STATES-1)
for ( i in 1:NUM_RATES ) {
    rate[i] := mu
}

Q := fnFreeBinary(rate, rescaled=false)

## references
[text](https://revbayes.github.io/tutorials/morph_ase/ase.html)
