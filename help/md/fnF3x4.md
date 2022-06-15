## name
fnF3x4

## title
The F3x4 codon frequency model

## description
This treats codon frequencies as a product of independent nucleotide frequencies.

Since stop codons are removed from the codon alphabet, frequencies are renormalized
so that the frequencies of non-stop codons sum to 1.0.

## details
## authors
## see_also
fnGY94, fnF1x4

## example
        kappa ~ dnLognormal(0,1)
        omega ~ dnUniform(0,1)
        pi1 ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
        pi2 ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
        pi3 ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
        Q := fnCodonGY94( kappa, omega, fnF3x4(pi1, pi2, pi3) )

## references
