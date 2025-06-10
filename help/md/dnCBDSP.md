## name
dnCBDSP
## title
Conditioned birth-death-shift process
## description
Simulates a tree under a birth-death process with shifts in birth and death
rates among lineages, conditioned on the assumption that no shifts take place
along extinct lineages.
## details
The initial birth and death rates can be specified with the `rootLambda` and
`rootMu` arguments, respectively. The rate at which speciation and extinction
rate shifts take place is specified by the `delta` argument, and the new
speciation or extinction rates are drawn from prior distributions specified in
the `lambda` and `mu` arguments. Similar to other birth-death processes in
RevBayes, `dnCBDSP` also takes arguments specifying the stopping `condition`
of the simulator (either survival or time) and the extant sampling probability
`rho`.

`dnCBDSP` is very similar to the model implemented in Bayesian Analysis
of Macroevolutionary Mixtures (BAMM; Rabosky 2014), particularly in making
a strong and potentially problematic assumption that all rate-shift events
have been observed (Moore et al. 2016) -- i.e., that no rate shifts are mapped
onto unobserved (extinct) branches. For an alternative birth-death-shift model
that relaxes this assumption, see `dnCDBDP` (Höhna et al. 2019), which employs
a finite number of rate categories instead of drawing rates directly from
a continuous distribution.
## authors
## see_also
dnCDBDP
mvBirthDeathEventContinuous
mvContinuousEventScale
mvEventTimeBeta
mvEventTimeSlide
## example
    # draw basic process parameters
    root_age ~ dnUniform(0, 2)
    root_lambda ~ dnUniform(0, 1)
    root_mu ~ dnUniform(0, 1)
    sampling_prob <- 1
    
    # simulate tree
    tree ~ dnCBDSP(rootAge    = root_age,
                   rootLambda = root_lambda,
                   rootMu     = root_mu,
                   delta      = 0.2,
                   rho        = sampling_prob,
                   condition  = "survival")
## references
- citation: Höhna S, Freyman WA, Nolen Z, Huelsenbeck JP, May MR, Moore BR (2019). A Bayesian approach for estimating branch-specific speciation and extinction rates. bioRxiv.
  doi: 10.1101/555805
  url: https://www.biorxiv.org/content/10.1101/555805v1.full
- citation: Moore BR, Höhna S, May MR, Rannala B, Huelsenbeck JP (2016). Critically evaluating the theory and performance of Bayesian analysis of macroevolutionary mixtures. Proceedings of the National Academy of Sciences of the USA, 113(34):9569-9574.
  doi: 10.1073/pnas.1518659113
  url: https://www.pnas.org/doi/full/10.1073/pnas.1518659113
- citation: Rabosky DL (2014). Automatic detection of key innovations, rate shifts, and diversity-dependence on phylogenetic trees. PLoS ONE, 9(2):e89543.
  doi: 10.1371/journal.pone.0089543
  url: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0089543
