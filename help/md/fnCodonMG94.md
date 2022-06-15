## name
fnCodonMG94

## title
The Muse-Gaut (1994) codon rate matrix

## description
The Muse-Gaut (1994) codon model.

A rate matrix on the 61 non-stop codons (in the standard genetic code).

Rates between codons with more than one nucleotide change are equal to 0.

In this model the rate Q(i,j) from i -> j is proportional to the frequency of
nucleotide in codon j that changed.  This differs from the Goldman-Yang (1994) model,
where Q(i,j) is proportional to the frequency of the entire codon j.

Unlike the Goldman-Yang (1994) model, the Muse-Gaut (1994) model does not allow all the codon
frequencies to vary independently.

## details
## authors
## see_also
fnCodonMG94, fnCodonMG94K

## example
        omega ~ dnUniform(0,1)
        pi ~ dnDirichlet( rep(2.0, 4) )
        Q1 := fnCodonMG94( omega, pi )

        Q2 := fndNdS( omega, fnX3( fnF81(pi) ) ) # MG94 = F81 + X3 + dNdS

## references
- citation: Muse, S. and B. Gaut (1994) A likelihood approach for comparing synonymous and nonsynonymous
      nucleotide substitution rates, with application to the chloroplast genome. Mol. Biol. Evol. (1994) 11 (5):715-724
  doi: https://doi.org/10.1093/oxfordjournals.molbev.a040152
   
