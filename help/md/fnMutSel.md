## name
fnMutSel

## title
Add mutation-selection balance to a rate matrix.

## description
Constructs a rate matrix from scaled selection coefficients w[i] and
mutation rate matrix mu(i,j).

A substitution from allele i -> j can be decomposed into
 (1) all individuals initially have state i
 (2) a single individual mutates from i -> j, at rate mu(i,j)
 (3) the allele j goes to fixation

Then the substitution rate Q is then given by
  Q(i,j) = mu(i,j) * Pr(j goes to fixation | i was fixed previously).

The probability of fixation is determined by scaled selection coefficients:
  w[i] = 2*N*s[i]
and the initial frequency 1/N of allele j.

## details
## authors
## see_also
fnCodonGY94, fnCodonMG94, fnX3, fndNdS

## example
        er ~ dnDirichlet( v(1,1,1,1,1,1) )
        nuc_pi ~ dnDirichlet( rep(2.0, 4) )
        w ~ dnIID(61, dnNormal(0,1))
        Q := fnMutSel(w, fnX3( fnGTR(er,nuc_pi) ) )   # GTR + X3 + MutSel

## references
- citation: Yang, Z. and R. Nielsen. Mutation-Selection Models of Codon Substitution and Their Use to Estimate
      Selective Strengths on Codon Usage.  Mol. Biol. Evol. (2008) 25(3):568--579
  doi: https://doi.org/10.1093/molbev/msm284
