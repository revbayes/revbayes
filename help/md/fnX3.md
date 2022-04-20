## name
fnX3

## title
Construct a codon rate matrix from a nucleotide rate matrix.

## description
Constructs a rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates of change from nucleotide i -> j at each codon position are given by the
nucleotide rate matrix.  The rate of 2 or 3 simultaneous changes is 0.

The X3 function can be used to construct other rate matrices in a modular fashion.
For example:
  (i)  MG94  = F81 + X3 + dNdS
  (ii) MG94K = HKY85 + X3 + dNdS

## details
## authors
## see_also
fnCodonGY94, fnCodonMG94K, fndNdS

## example

        kappa ~ dnLognormal(0,1)
        omega ~ dnUniform(0,1)
        nuc_pi ~ dnDirichlet( rep(2.0, 4) )
        Q1 := fnCodonMG94K( kappa, omega, nuc_pi )
        # This is the same.
        Q2 := fndNdS(fnX3(fnHKY(kappa, nuc_pi)), omega)          # HKY + X3 + dNdS, or HKY*3 + dNdS

        er ~ dnDirichlet( v(1,1,1,1,1,1) )
        Q3 := fnX3(fnGTR(er, nuc_pi))      # GTR + X3, or GTR*3

## references
- citation: Redelings, BD (2021). BAli-Phy version 3: Model-based co-estimation of Alignment
       and Phylogeny.  Bioinformatics (2021) 37(10):3032–3034.
  doi: https://doi.org/10.1093/bioinformatics/btab129
