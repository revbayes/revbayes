## name
fnCodonMG94K

## title
The Muse-Gaut (1994) codon rate matrix + K.

## description
The Muse-Gaut (1994) codon model, extended with a transition/transversion rate ratio.

A rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates between codons with more than one nucleotide change are equal to 0.

In this model the rate Q(i,j) from i -> j is proportional to the frequency of
nucleotide in codon j that changed.  This differs from the Goldman-Yang (1994) model,
where Q(i,j) is proportional to the frequency of the entire codon j.

This version is an extension of the fnCodonMG94 model to add a transition/transversion
rate ratio.  This makes it more comparable to the Goldman-Yang (1994) model.

Unlike the Goldman-Yang (1994) model, the Muse-Gaut (1994) model does not allow all the codon
frequencies to vary independently.

## details
## authors
## see_also
fnCodonGY94, fnCodonMG94K

## example
        kappa ~ dnLognormal(0,1)
        omega ~ dnUniform(0,1)
        pi ~ dnDirichlet( rep(2.0, 4) )
        Q1 := fnCodonMG94K( kappa, omega, pi )

        Q2 := fndNdS( omega, fnX3( fnHKY( kappa, pi) ) ) # MG94K = HKY + X3 + dNdS

## references
- citation: Muse, S. and B. Gaut (1994) A likelihood approach for comparing synonymous and nonsynonymous
      nucleotide substitution rates, with application to the chloroplast genome. Mol. Biol. Evol. (1994) 11 (5):715-724
  doi: https://doi.org/10.1093/oxfordjournals.molbev.a040152
   
