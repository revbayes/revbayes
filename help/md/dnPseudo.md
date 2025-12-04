## name
dnPseudo

## title
Pseudodistribution

## description
Fake distribution for adding likelihood evidence

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
hypothetical observation F such that

    Pr(F|t, a) = 1 if t > a
               = 0 otherwise

For example, F might represent a fossil with a stratigraphic range.
Conceptually, there exists some (unknown and irrelevant) distribution D
for F such that

    F ~ D(t,a)

and we only want to specify its likelihood function Pr(F|t, a).

In this case, we say that F is "pseudodata"; it is sometimes called
"virtual evidence" in the literature.  We write this as

    F ~ dnPseudo(t, pdAbove(a))
    F.clamp( pseudoObservation() )

Here pdAbove(a) specifies the likelihood function.  We call it the
pseudodata likelihood. It is NOT a probability distribution, because it
is not normalized to sum to 1.  The likelihood at x for a pseudodata
likelihood pd can be obtained with pd.lnLikelihood(x).

Pseudodata likelihoods can be combined in several ways:

  * pd1 && pd2: multiplies the likelihoods.
        For example pdAbove(a) && pdBelow(b) = pdBetween(a,b)

  * pd1 || pd2: adds the likelihoods.

  * pd + d: shifts the likelihood curve right by d.  
        For example, pdLogNormal(0,1) + 10 requires that x > 10.

  * pd - d: shifts the likelihood curve left by d.

  * d - pd: computes pd(d-x).

## authors
Benjamin Redelings

## see_also
pdWeightReal
pdBetween
pdAbove
pdBelow
pdNormal
pdExponential
pdLogNormal

## example
    # A taxon age range
    f ~ dnPseudo(t[i], pdBetween(a[i], b[i]))
    f.clamp( pseudoObservation() )

    # A fossil node calibration
    f ~ dnPseudo(t[i], pdLogNormal(lmu, lsigma) + shift)  
    f.clamp()

    # Observing that taxon 1 is younger than taxon 2
    f ~ dnPseudo(taxon1_age, pdAbove(taxon2_age))
    f.clamp( pseudoObservation() )

    # Checking likelihood values
    pdBetween(0, 1).lnLikelihood(0.5)  # yields 0
    pdBetween(0, 1).lnLikelihood(1.5)  # yields -inf

    # A Mixture likelihood
    (pdWeightReal(p) && pd1) || (pdWeightReal(1-p) && pd2)

    # A soft constraint    
    pdWeightReal(1/1000) || pd
