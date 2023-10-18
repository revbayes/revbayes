## name
fnF2x4

## title
The F2x4 doublet frequency model

## description
This treats doublet frequencies as a product of independent nucleotide frequencies.

## details
## authors
## see_also
fnX2

## example
        # An RNA stem model with independent base frequencies (from fnF2x4),
        # and simultaneous 2-nucleotide changes allows.
        nuc_pi ~ dnDirichlet( v(2.0, 2.0, 2.0, 2.0) )
        rna_stem_er ~ dnDirichlet( rep(1.0, 16*15/2) )
        rna_stem_pi := fnF2x4(nuc_pi, nuc_pi)
        Q := fnGTR(rna_stem_er, rna_stem_pi)

## references
