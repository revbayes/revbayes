## name
fnF1x4

## title
The F1x4 codon frequency model

## description
This treats codon frequencies as a product of independent nucleotide frequencies.

Since stop codons are removed from the codon alphabet, frequencies are renormalized
so that the frequencies of non-stop codons sum to 1.0.

## details
## authors
## see_also
fnGY94, fnF3x4

## example
        kappa ~ dnLognormal(0,1)
        omega ~ dnUniform(0,1)
        pi ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
        Q := fnCodonGY94( kappa, omega, fnF1x4(pi) )

## references
