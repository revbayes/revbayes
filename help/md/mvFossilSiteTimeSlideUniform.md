## name
mvFossilSiteSlideUniform
## title
Move to simultaneously update the ages of site-linked fossils
## description
Proposes simultaneous additive updates to the ages of fossil taxa that come
from the same site and are therefore known to be contemporaneous, even if their
exact age remains uncertain (King & Rücklin 2020).
## details
The move takes a vector of `Taxon` objects, representing the fossils whose ages
are to be treated as linked. The minimum and maximum ages of the site that the
fossils come from can either be supplied to the move via the `min` and `max`
arguments, or taken directly from the `Taxon` objects. Internally, the minimum
and maximum values are recomputed to account for the current tree topology,
which may impose additional constraints beyond those provided by the age range
of the site (since a node cannot be older than its parent, and a sampled
ancestor cannot be younger than its oldest descendant). Like other slide moves,
`mvFossilSiteSlideUniform` then draws a uniformly distributed random number and
adds the draw to the current value. The width of the uniform distribution
determines the "boldness" of the proposal; it can be specified by the `delta`
argument and, if `tune=TRUE`, automatically adjusted to ensure that the
acceptance rate of the move reaches `tuneTarget`. The minimum and maximum ages
are used as reflection boundaries if the addition of the uniform draw to the
current age produces a value that is too low or too high.

Note that there are cases in which the ages of site-linked taxa may (seemingly)
diverge even if they are updated using `mvFossilSiteSlideUniform`. This may
occur when the ages of the taxa in question are also adjusted using moves that
are not site-aware (`mvFossilTipTimeSlideUniform`, `mvFossilTipTimeUniform`).
These moves can be used without specifying the `tip` argument, in which case
they randomly pick a fossil every time they are performed, and may consequently
update the age of one site-linked fossil but not of the others. Similarly, if
an analysis also employs the `mvCollapseExpandFossilBranch` move, which turns
tips into sampled ancestors and _vice versa_, the resulting summary trees may
show the site-linked taxa to have different ages even if their ages are in fact
exactly identical in every single individual MCMC sample. This is because the
corresponding point estimates are averaged over different subsets of samples,
determined by the status of a given fossil as a tip or a sampled ancestor. The
behavior can be suppressed by specifying `conditionalAges=TRUE` in the tree
summary functions (`mapTree()`, `mccTree()`), which ensures that the tip age
point estimates are averaged across only those samples whose entire topology
equals the summary tree.
## authors
Mario Schädel
## see_also
mvCollapseExpandFossilBranch
mvFossilTipTimeSlideUniform
## example
    # Assume that the min/max ages of each taxon have been read from a TSV file
    # Extract all fossil taxa in the tree
    fossils = fbd_tree.getFossils()
    
    # Specify a uniform prior spanning the age uncertainty range of each fossil
    for(i in 1:fossils.size())
    {
        t[i] := tmrca(fbd_tree, clade(fossils[i]))
        
        a[i] = fossils[i].getMinAge()
        b[i] = fossils[i].getMaxAge()
        
        F[i] ~ dnUniform(t[i] - b[i], t[i] - a[i])
        F[i].clamp( 0 )
    }
    
    # Create a vector including all fossils known from the same site
    sitelinked_fossils[1] = taxon("fossil_X")
    sitelinked_fossils[2] = taxon("fossil_Y")
    sitelinked_fossils[3] = taxon("fossil_Z")
    
    # Add a move to update the ages of these fossils simultaneously
    moves.append( mvFossilSiteSlideUniform(fbd_tree, origin_time, taxa=sitelinked_fossils,
                                           weight=10.0) )
                                           
    # The ages of the remaining fossils should be updated independently
    nonsitelinked_fossils = fossils
    nonsitelinked_fossils.erase( sitelinked_fossils )
    
    for(i in 1:nonsitelinked_fossils.size())
    {
        moves.append( mvFossilTimeSlideUniform(fbd_tree, origin_time,
                                               tip=nonsitelinked_fossils[i], weight=5.0) )
    }
## references
- citation: King B, Rücklin M (2020). Tip dating with fossil sites and stratigraphic sequences. PeerJ, 8:e9368.
  doi: 10.7717/peerj.9368
  url: https://peerj.com/articles/9368/
