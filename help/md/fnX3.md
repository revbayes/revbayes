## name
fnX3

## title
Construct a codon rate matrix from a nucleotide rate matrix.

## description
Constructs a rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates of change from nucleotide i -> j at each codon position are given by the
nucleotide rate matrix.  The rate of 2 or 3 simultaneous changes is 0.

## details
## authors
## see_also
fnCodonGY94, fnCodonMG94K

## example

        kappa ~ dnLognormal(0,1)
        nuc_pi ~ dnDirichlet( rep(2.0, 4) )
        Q1 := fnCodonMG94K( kappa, 1.0, nuc_pi )
        # This is the same.
        Q2 := fnX3(fnHKY(k,nuc_pi))   # HKY + X3, or HKY*3

        er ~ dnDirichlet( v(1,1,1,1,1,1) )
        Q3 := fnX3(fnGTR(er,pi))      # GTR + X3, or GTR*3

## references
   
