## name
dnPseudo

## title
Pseudodistribution

## description
Virtual distribution for adding likelihood evidence

## details
Researchers often want to add information about a parameter that
already has a prior. For example, we might want to enforce that a
taxon age t is in the range [a, infinity) or [a,b].  It might seem
natural to specify a "second prior", but this does not work:

  * A uniform distribution on [a,infinity) cannot be normalized
    to sum to 1.
  
  * A uniform distribution on [a,b] can be normalized, but 
    normalizing it leads to incorrect results if a and b are
    themselves random.

The correct approach is to treat this extra information as likelihood
information.  If we know that t > a, then we can introduce a 
hypothetical observation F=f such that

    Pr(F=f|t, a) = 1 if t > a
                 = 0 otherwise

For example, F could be the stratigraphic depth for a fossil.
Conceptually, there exists some distribution D for F such that

    F ~ D(t,a)

Note that the stratigraphic depth is considered to be a random
variable that depends on the age t, and not the other way around.

Now, we may want to avoid explicitly modeling the stratigraphic
observation and simply specify the likelihood function Pr(F=f|t, a)
directly.  In this case the observed depth f and the distribution
D are unspecified.

In this case, we say that F=f is "pseudodata" or "virtual evidence".
We write this as

    F ~ dnPseudo( pdAbove(t,a) )
    F.clamp( pseudoObservation() )

Here dnPseudo takes the place of the unspecified distribution D,
and pseudoObservation() takes the place of the unspecified
stratigraphic depth.  The only thing that is specified is the
likelihood function pdAbove(t,a).  We call it a pseudodata likelihood
because the true observed value and model are hidden.

It is important to note that his is NOT a probability distribution,
because it is not normalized to sum to 1.

Pseudodata likelihoods can be combined by the && operator.  This
operator multiplies the likelihoods. For example 

    pdAbove(t,a) && pdBelow(t,b)

Is the same as 

    pdBetween(t,a,b)

## authors
Benjamin Redelings

## see_also
pdAbove
pdBelow
pdBetween
pdNormal
pdExponential
pdLogNormal
pdConstant
pdRequire

## example
    # A taxon age range
    f ~ dnPseudo( pdBetween(t[i], a[i], b[i]))
    f.clamp( pseudoObservation() )

    # A fossil node calibration
    f ~ dnPseudo( pdLogNormal(t[i], lmu, lsigma, shift) )
    f.clamp()

    # Observing that taxon 1 is younger than taxon 2
    f ~ dnPseudo( pdBelow(taxon1_age, taxon2_age) )
    f.clamp( pseudoObservation() )

    # Checking likelihood values
    pdBetween(0.5, 0, 1)  # yields PDL0.0
    pdBetween(1.5, 0, 1)  # yields PDL-inf
