## name
mvFossilTipTimeUniform
## title
Move to uniformly draw fossil tip ages
## description
This moves either takes a specific fossil, or randomly picks a fossil, and then draws the new ages randomly between the maximum and minimum ages.

## details
The maximum ages is computed either by its parents or the maximum age in the uncertainty of the fossil, which can be provided to the move or is taken from the taxon object.
The minimum ages is computed either by its oldest descendant (for sampled ancestors) or the minimum age in the uncertainty of the fossil, which can be provided to the move or is taken from the taxon object.

## authors
Sebastian Hoehna

## see_also
mvFossilTipTimeSlideUniform

## example

# Use a for loop to create a uniform distribution on the occurrence time for each fossil #
# The boundaries of the uniform distribution are specified in the tsv file #
fossils = fbd_tree.getFossils()
for(i in 1:fossils.size())
{
    t[i] := tmrca(fbd_tree, clade(fossils[i]))

    a[i] = fossils[i].getMinAge()
    b[i] = fossils[i].getMaxAge()

    F[i] ~ dnUniform(t[i] - b[i], t[i] - a[i])
    F[i].clamp( 0 )
    moves.append( mvFossilTipTimeUniform(fbd_tree, origin_time, min=a[i], max=b[i], tip=fossils[i], weight=5.0) )
    moves.append( mvFossilTipTimeSlideUniform(fbd_tree, origin_time, min=a[i], max=b[i], tip=fossils[i], weight=5.0) )
}


## references
