## name
fnCodonGY94

## title
The Goldman-Yang (1994) codon rate matrix

## description
The Goldman-Yang (1994) codon model.

A rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates between codons with more than one nucleotide change are equal to 0.

In this model the rate Q(i,j) from i -> j is proportional to the frequency of codon j.
This means that the rate of change between low-frequency codons is lower than the rate
between high-frequency codons, even when the nucleotide change involved is the same.
In other words, the rate of change from nucleotide n1 -> n2 depends on its neighboring
nucleotides.  This differs from the Muse-Gaut (1994) model, and is perhaps less realistic.

Unlike the Muse-Gaut (1994) model, the Goldman-Yang (1994) model can allow all the codon
frequencies to vary independently.

## details
## authors
## see_also
fnF1x4, fnF3x4, fnCodonMG94, fnCodonMG94K

## example
        kappa ~ dnLognormal(0,1)
        omega ~ dnUniform(0,1)
        pi61 ~ dnDirichlet( rep(2.0, 61) )
        Q1 := fnCodonGY94( kappa, omega, pi61 )

        pi1 ~ dnDirichlet( rep(2.0, 4) )
        Q2 := fnCodonGY94( kappa, omega, fnF1x4(pi1) )

        pi2 ~ dnDirichlet( rep(2.0, 4) )
        pi3 ~ dnDirichlet( rep(2.0, 4) )
        Q3 := fnCodonGY94( kappa, omega, fnF3x4(pi1, pi2, pi3) )

## references
- citation: Goldman, N. and Z. Yang (1994). A codon-based model of nucleotide substitution for protein-coding DNA
      sequences. Mol. Biol. Evol. (1994) 11 (5):725-736
  doi: https://doi.org/10.1093/oxfordjournals.molbev.a040153
