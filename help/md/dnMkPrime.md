## name
MkPrime
## title
Mk-prime model: phylogenetic CTMC with inference of unobserved character states
## description
Computes the probability of observed discrete character data under the Mk'
(Mk-prime) model, which extends the standard Mk model by inferring the number
of unobserved character states `u` beyond those observed in the data.

Characters must be grouped into a partition sharing the same observed state
count `kObs`. The total number of states is `k = kObs + u`. The substitution
process is a JC(k) continuous-time Markov chain; the relabelling correction
ensures the likelihood accounts for the unidentifiability of state labels.

The prior on `u` is specified internally. Common choices are geometric
(default), log-series, or power-law.
## details
The log-likelihood computed by `dnMkPrime` is:

  log L = log L_CTMC(kObs + u) + N * log C(kObs + u, kObs) + log p(u)

where:
- `log L_CTMC(k)` is the standard Mk likelihood under a JC(k) matrix,
  summed over patterns weighted by pattern counts and conditioned on the
  ascertainment bias correction specified by `coding`,
- `N` is the number of characters in the partition,
- `log C(k, kObs)` is the relabelling correction
  (log kObs! - log k! + kObs*log(k)),
- `log p(u)` is the prior on u (geometric, logseries, or powerlaw).

At `u = 0`, the relabelling correction is 0 and the likelihood reduces to the
standard conditioned Mk likelihood, equal to `dnPhyloCTMC` with `coding="variable"`.

`u` must be declared as a stochastic node in the Rev script (e.g., via
`u ~ dnGeometric(0.5)`) so that standard MCMC moves can sample it. The
external prior on `u` should be flat (e.g., `dnUniform(0, kMax - kObs)`) since
`dnMkPrime` applies the internal prior internally; using an informative external
prior will double-count.

Note: on short trees (total branch length << 1), the relabelling correction
dominates and the posterior on `u` will concentrate at 0 regardless of `kObs`.
Inference of `u > 0` is most reliable when total tree length is approximately
1-2 substitutions per site.
## authors
## see_also
dnPhyloCTMC
## example
    # Read morphological character data
    morph <- readDiscreteCharacterData("morpho.nex")

    # Define taxa and draw a tree
    taxa <- morph.taxa()
    tree ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExp(10))

    # Partition characters by observed state count
    morph_bin <- morph
    morph_bin.setNumStatesPartition(2)
    morph_bin.includeCharacter(morph_bin.getIncludedCharacterIndices())

    # Declare unobserved state count (stochastic, with flat external prior)
    u ~ dnUniformNatural(0, 6)
    moves.append(mvRandomNaturalWalk(u, weight=1))

    # Clamp data to the Mk-prime distribution
    chars ~ dnMkPrime(
        tree    = tree,
        kObs    = 2,
        u       = u,
        prior   = "geometric",
        priorParam = 0.4,
        coding  = "variable"
    )
    chars.clamp(morph_bin)

