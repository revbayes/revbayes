## name
fnFoldSFS
## title
Fold a site frequency spectrum (SFS)
## description
Converts an unfolded site frequency spectrum with (n + 1) elements for
n individuals to a folded site frequency spectrum with (n/2 + 1) elements.
## details
An unfolded site frequency spectrum for n individuals can be represented as
a vector in which the first element gives the number of monomorphic sites at
which all sampled individuals carry the ancestral state, the next (n - 1)
elements give the number of polymorphic sites at which the derived allele is
present in 1, ..., n - 1 individuals, and the last element gives the number
of monomorphic sites at which all the individuals carry the derived state.
This gives a total of (n + 1) elements for n individuals.

If we do not know which allele is ancestral and which is derived, we can only
distinguish between minor (rarer) and major (more common) alleles. In this
case, we "fold" the spectrum to obtain the counts of minor alleles by grouping
together counts k and (n + 2 - k): e.g., for 10 individuals, we would collapse
together entries 1 and 11 (the monomorphic sites), 2 and 10, 3 and 9, 4 and 8,
and 5 and 7. The entry k = n/2 + 1 (i.e., the 6th element in this example),
which gives the number of sites at which the allele frequencies are exactly
50/50, is unique and not grouped together with any other entry.
## authors
Sebastian Höhna
David Černý
## see_also
## example
    # let's assume we have some SFS "observed"
    obs_sfs = [ 305082, 44248, 32223, 28733, 28220, 26205, 27477, 26618, 27533,
                26945, 28736, 28671, 31277, 31250, 34352, 34859, 38331, 40005,
                45666, 48986, 65829, 64363, 70895, 74114, 82705, 88226, 102194,
                114566, 130176, 143775, 169216, 191624, 230016, 276489, 333069,
                394810, 501961, 653809, 890077, 1349350, 50296796 ]
                
    folded_sfs = fnFoldSFS( obs_sfs )
## references
