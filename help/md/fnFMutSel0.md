## name
fnFMutSel0

## title
The FMutSel0 model

## description
Constructs a rate matrix from 61 scaled selection coefficients w[i] and
a 4x4 nucleotide mutation rate matrix mu(i,j).  In the original paper
the nucleotide mutation rate matrix is a GTR rate matrix.

The FMutSel0 model is a restriction of the FMutSel model that constrains
all codons for the same amino acid to have the same scaled selection
coefficient.

The function fnMutSelAA differs from fnFMutSel0 by taking a codon mutation
rate matrix.

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
fnCodonGY94, fnCodonMG94, fnFMutSel0, fnMutSel

## example
        er ~ dnDirichlet( v(1,1,1,1,1,1) )
        nuc_pi ~ dnDirichlet( rep(2.0, 4) )
        F ~ dnIID(20, dnNormal(0,1))
        omega ~ dnUniform(0,1)
        # The FMutSel0 model from Yang and Nielsen (2008)        
        Q1 := fnFMutSel0(fnGTR(er, nuc_pi), F, omega)

        # The same -- fMutSel0 = GTR(er,nuc_pi) + X3 + MutSel(F) + dNdS(omega)
        Q2 := fndNdS( fnMutSelAA( fnX3( fnGTR(er, nuc_pi)), F), omega)

## references
- citation: Yang, Z. and R. Nielsen. Mutation-Selection Models of Codon Substitution and Their Use to Estimate
      Selective Strengths on Codon Usage.  Mol. Biol. Evol. (2008) 25(3):568--579
  doi: https://doi.org/10.1093/molbev/msm284
