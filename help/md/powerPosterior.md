## name
powerPosterior
## title
Power posterior analysis
## description
Samples from a series of "power posterior" distributions with the likelihood
term raised to a power between 0 and 1. Such distributions are often used to
estimate marginal likelihoods for model selection or hypothesis testing via
Bayes factors (Gelman & Meng 1998; Friel & Pettitt 2008).
## details
A power posterior analysis samples from a series of importance distributions of
the form:
    f_beta(theta | Y) = f(theta) x f(Y | theta)^beta
    
where theta jointly denotes the parameters of interest and Y denotes data, so
f(theta) denotes the prior and f(Y | theta) denotes the likelihood. Since beta
ranges from 0 to 1, the power posterior distributions form "stepping stones"
along the path between the prior (beta = 0) and the posterior (beta = 1). For
this reason, individual power posteriors are also referred to as stones, and
the two main techniques that employ them to estimate marginal likelihoods are
known as path sampling (Gelman & Meng 1998; Lartillot & Philippe 2006) and
stepping-stone sampling (Fan et al. 2011; Xie et al. 2011).

The user can either supply their own vector of beta powers using the `powers`
argument, or just a number of power posteriors to sample from using the `cats`
argument. In the latter case, if `cats=K`, the powers are calculated following
Xie et al. (2011) as:
    beta_i = [i / (K - 1)]^(1 / alpha) for i in K - 1, ..., 0
    
so that they correspond to evenly spaced quantiles of the Beta(alpha, 1)
distribution, where the shape parameter `alpha` is a user-specified argument
set to 0.2 by default. The samples from each distribution are recorded in a
file with a base name specified by the `filename` argument and an automatically
appended suffix equal to (K - i). The output file with the suffix "_stone_1" 
will therefore contain samples from the actual posterior (beta = 1), while the
file with the suffix "_stone_<K>" will contain samples from the prior. For 
reasons of numerical stability, the beta of this final stone is not set exactly
to 0 but to an extremely small positive number (~ 1.2e-302).

The RevBayes implementation of power posterior analysis is fully parallelized
(Hoehna et al. 2021). In general, after a common pre-burnin stage that should
allow the sampler to converge to the posterior (see the `.burnin()` method),
the `.run()` method will distribute the K stones among M available CPUs in such
a way that each CPU handles floor(K/M) or ceiling(K/M) consecutive powers. This
has the advantage of allowing the last sample for one power to be used as the
starting state for the subsequent power. However, to account for the transition
from one power to the next, each power posterior should still include a small
burnin fraction (set to 0.25 by default and specified by the `burninFraction`
argument to the `.run()` method). For example, with K = 50, M = 8, and a single
CPU used for each likelihood computation (`procPerLikelihood=1`, by default),
the individual power posteriors will be distributed among the CPUs as follows:

    -----------------------------------------------------------------
                        Filename suffix (= K - i)
    -----------------------------------------------------------------
    | CPU 1 | CPU 2 | CPU 3 | CPU 4 | CPU 5 | CPU 6 | CPU 7 | CPU 8 |
    |-------|-------|-------|-------|-------|-------|-------|-------|
    |   1   |   7   |   13  |   19  |   26  |   32  |   38  |   44  |
    |   2   |   8   |   14  |   20  |   27  |   33  |   39  |   45  |
    |  ...  |  ...  |  ...  |  ...  |  ...  |  ...  |  ...  |  ...  |
    |   6   |   12  |   18  |   24  |   31  |   37  |   43  |   49  |
    |       |       |       |   25  |       |       |       |   50  |

More flexible strategies are enabled by the `.runOneStone()` method, which
allows the user to execute a separate analysis for each power posterior.
## authors
Sebastian Hoehna
Michael Landis
John Huelsenbeck
## see_also
pathSampler
steppingStoneSampler
## example
    # Create a simple model (unclamped)
    a ~ dnExponential(1)
    mymodel = model(a)
    
    # Create a move vector and a monitor vector
    moves[1] = mvScale(a, lambda = 1.0, weight = 1.0)
    monitors[1] = mnFile(a, filename = "output/out.log")
    
    # Create an analysis object to sample from 16 distributions
    pow_p = powerPosterior(mymodel, monitors, moves, "output/out.pp", cats=16, sampleFreq=1)
    
    # Execute a single power posterior ("stone"), or the entire analysis
    pow_p.burnin(generations=100, tuningInterval=50)                # pre-burnin
    pow_p.runOneStone(index=1, generations=20, burninFraction=0.1)  # run just the posterior
    pow_p.runOneStone(index=16, generations=20, burninFraction=0.1) # run just the prior
    pow_p.run(generations=20, burninFraction=0.1)                   # run all stones
    
    # Compute the marginal likelihood using the stepping-stone sampler
    ss = steppingStoneSampler(file="output/out.pp", powerColumnName="power",
                              likelihoodColumnName="likelihood")
    ss.marginal()
## references
- citation: Fan Y, Wu R, Chen M-H, Kuo L, Lewis PO (2011). Choosing among partition models in Bayesian phylogenetics. Molecular Biology and Evolution, 28(1):523-532.
  doi: 10.1093/molbev/msq224
  url: https://academic.oup.com/mbe/article/28/1/523/983866
- citation: Friel N, Pettitt AN (2008). Marginal likelihood estimation via power posteriors. Journal of the Royal Statistical Society Series B: Statistical Methodology, 70(3):589-607.
  doi: 10.1111/j.1467-9868.2007.00650.x
  url: https://academic.oup.com/jrsssb/article-abstract/70/3/589/7109555
- citation: Gelman A, Meng X-L (1998). Simulating normalizing constants: from importance sampling to bridge sampling to path sampling. Statistical Science, 13(2):163-185.
  doi: 10.1214/ss/1028905934
  url: https://www.jstor.org/stable/2676756
- citation: Hoehna S, Landis MJ, Huelsenbeck JP (2021). Parallel power posterior analyses for fast computation of marginal likelihoods in phylogenetics. PeerJ, 9:e12438.
  doi: 10.7717/peerj.12438
  url: https://peerj.com/articles/12438/
- citation: Lartillot N, Philippe H (2006). Computing Bayes factors using thermodynamic integration. Systematic Biology, 55(2):195–207.
  doi: 10.1080/10635150500433722
  url: https://academic.oup.com/sysbio/article-abstract/55/2/195/1620800
- citation: Xie W, Lewis PO, Fan Y, Kuo L, Chen M-H (2010). Improving marginal likelihood estimation for Bayesian phylogenetic model selection. Systematic Biology, 60(2):150-160.
  doi: 10.1093/sysbio/syq085
  url: https://academic.oup.com/sysbio/article-abstract/60/2/150/2461669
