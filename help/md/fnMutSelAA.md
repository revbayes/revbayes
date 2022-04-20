## name
fnMutSelAA

## title
Add mutation-selection balance to a rate matrix -- fitnesses on amino acids

## description
Constructs a rate matrix from scaled selection coefficients w[i] and
mutation rate matrix mu(i,j).

fnMutSelAA takes 20 scaled selection coefficients, one for each amino acid.
This differs from fnMutSel, which takes 61 scaled selection coefficients,
one for each codon.  fnMutSelAA assumes that codons for the same amino acid
have the same fitness.

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
fnCodonGY94, fnCodonMG94, fnX3, fndNdS, fnMutSel

## example
        er ~ dnDirichlet( v(1,1,1,1,1,1) )
        nuc_pi ~ dnDirichlet( rep(2.0, 4) )
        F ~ dnIID(20, dnNormal(0,1))
        Q := fnMutSelAA(fnX3(fnGTR(er, nuc_pi)), F)

## references
- citation: Yang, Z. and R. Nielsen. Mutation-Selection Models of Codon
      Substitution and Their Use to Estimate Selective Strengths on Codon
      Usage.  Mol. Biol. Evol. (2008) 25(3):568--579
  doi: https://doi.org/10.1093/molbev/msm284
