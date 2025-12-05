## name
dnFastBirthDeathShift
## title
The relaxed birth-death-shift process
## description
Simulates a tree under the four-parameter birth-death-shift process 
## details
This distribution represents the birth-death-shift process, which has four
possible events:

- a speciation event (lambda)
- an extinction event (mu)
- a speciation rate shift event (alpha)
- an extinction rate shift event (beta)

The possible speciation rates (lambda) and extinction rates (mu) are
assumed to be distributed according to discretized base distributions.
Currently, the log-normal base distribution is hard coded. The log-normal
base distributions have location parameters (speciationScale, extinctionScale) 
and spread parameters (speciationSD, extinctionSD), which can be fixed or
estimated using prior distributions.

In comparison with dnCDBDP, the dnFastBirthDeathShift distribution is 
simplified and more optimized, allowing it to be significantly faster than 
the implementation in dnCDBDP (some tests indicate 60-70x performance).

The model resembles the Birth-Death-Shift model of Hoehna et al. (2019), 
but there are two differences: 

1. speciation rate shift events and extinction rate shift events are allowed 
   to occur at different rates (i.e., separate parameters alpha and beta 
   instead of a single eta).
2. rate shift events that lead to a simultaneous change in the speciation 
   and extinction rate are not allowed.

To estimate the model parameters using data, use the .clamp() method to fix
a phylogeny. The clamped phylogeny must be in units of time, only have extant 
tips, have non-negative branch lengths and be bifurcating. Non-branching 
internal nodes (i.e., sampled ancestors) are not allowed.


## authors
Bjorn Tore Kopperud
## see_also
dnCBDSP
## example
    # specify the tree height
    root = 8.0

    # Specify the magnitude of the base distributions
    speciation_mean <- 0.2
    extinction_mean <- 0.15

    # The spread of the base distributions
    speciation_sd <- 0.6
    extinction_sd <- 0.6

    # number of discretizations in the base distributions
    NUM_RATE_CLASSES <- 6

    # the parameters for the number of rate shift events
    speciation_shift_rate <- 0.02
    extinction_shift_rate <- 0.015


    # construct a variable for the tree drawn from a birth death shift process
    timetree ~ dnFastBirthDeathShift( rootAge            = root,
                                    speciationScale      = speciation_mean,
                                    extinctionScale      = extinction_mean,
                                    speciationSD         = speciation_sd,
                                    extinctionSD         = extinction_sd,
                                    alpha                = speciation_shift_rate,
                                    beta                 = extinction_shift_rate,
                                    numSpeciationClasses = NUM_RATE_CLASSES,
                                    numExtinctionClasses = NUM_RATE_CLASSES,
                                    rho                  = rho,
                                    condition            = "time")


## references
- citation: Barido-Sottani J, Vaughan TG, Stadler T (2020). A multitype birth–death model for Bayesian inference of lineage-specific birth and death rates. Systematic Biology, 69(5):973–986.
  doi: 10.1093/sysbio/syaa016
  url: https://academic.oup.com/sysbio/article/69/5/973/5762626
- citation: Höhna S, Freyman WA, Nolen Z, Huelsenbeck JP, May MR, Moore BR (2019). A Bayesian approach for estimating branch-specific speciation and extinction rates. bioRxiv.
  doi: 10.1101/555805
  url: https://www.biorxiv.org/content/10.1101/555805v1.full
