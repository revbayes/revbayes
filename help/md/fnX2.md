## name
fnX2

## title
Construct a doublet (16x16) rate matrix from a nucleotide rate matrix.

## description
Constructs a double rate matrix on the 16 nucleotide pairs.

Rates of change from nucleotide i -> j at each doublet position are given by the
nucleotide rate matrix.  The rate of 2 simultaneous changes is 0.

The X3 function can be used to constructor rate matrices on doublets in a
modular fashion.

## details
## authors
## see_also
fnX3

## example

        kappa ~ dnLognormal(0,1)
        nuc_pi ~ dnDirichlet( rep(2.0, 4) )
        # Mutation rate matrix on RNA stems
        Q1 := fnX2( fnHKY(kappa, nuc_pi) )
        F ~ dnIID(16, dnNormal(0,1))
        # Add selection to the rate matrix
        Q2 := fnMutSel(Q1, F)

## references
