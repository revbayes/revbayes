## name
fnMutSel

## title
Add mutation-selection balance to a rate matrix.

## description
Constructs a rate matrix from scaled selection coefficients w[i] and
mutation rate matrix mu(i,j).

fnMutSel takes 61 scaled selection coefficients, one for each codon.
This differs from fnMutSelAA, which takes 20 scaled selection coefficients,
one for each amino acid.

A substitution from allele i -> j can be decomposed into
 (1) all individuals initially have state i
 (2) a single individual mutates from i -> j, at rate mu(i,j)
 (3) the allele j goes to fixation

Then the substitution rate Q is then given by
  Q(i,j) = mu(i,j) * Pr(j goes to fixation | i was fixed previously).

The probability of fixation is determined by scaled selection coefficients:
  F[i] = 2*N*s[i]
and the initial frequency 1/N of allele j.

## details
## authors
## see_also
fnCodonGY94, fnCodonMG94, fnMutSelAA, fnFMutSel, fndNdS

## example
        er ~ dnDirichlet( v(1,1,1,1,1,1) )
        nuc_pi ~ dnDirichlet( rep(2.0, 4) )
        F ~ dnIID(61, dnNormal(0,1))
        Q := fnMutSel(fnX3(fnGTR(er, nuc_pi) ), F)       # GTR + X3 + MutSel

        # A mutation-selection balance model on RNA, with GTR mutation.
        F2 ~ dnIID(16, dnNormal(0,1))
        Q2 := fnMutSel(fnX2(fnGTR(er,nuc_pi) ), F2)      # GTR + X2 + MutSel

## references
- citation: Yang, Z. and R. Nielsen. Mutation-Selection Models of Codon Substitution and Their Use to Estimate
      Selective Strengths on Codon Usage.  Mol. Biol. Evol. (2008) 25(3):568--579
  doi: https://doi.org/10.1093/molbev/msm284
