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
Description

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

Implementation

A power posterior analysis consists of up to four distinct stages:

    1. Global pre-burnin
    2. Stone-specific pre-burnin
    3. Stone-specific burnin
    4. Stone-specific sampling phase

The global pre-burnin (1) is executed using the `.burnin()` method and is meant
to allow the sampler to converge to the posterior. No samples are collected, 
and a separate generation counter is used, much like with `mcmc.burnin()`. This
stage is executed using a single worker, which corresponds either to a single
process, or to a group of `procPerLikelihood` processes if this argument
(controlling the number of processes used for each likelihood computation) is
set to a value other than 1 (default).

In contrast, the remaining three stages are executed using the `.run()` method
and fully parallelized (see Höhna et al. 2021, though note that the current
implementation diverges somewhat from that described therein). The K stones are
distributed among M available workers in such a way that each worker handles
floor(K/M) or ceiling(K/M) consecutive powers. This allows the last sample for
one power to be used as the starting state for the next power. The assignment
of stones to available workers is summarized in a table printed at the start of
every analysis. For example, with K = 50, M = 8, and `procPerLikelihood=1`, the
individual power posteriors will be distributed among processes as follows:

    Process | Step 1 | Step 2 | Step 3 | Step 4 | Step 5 | Step 6 | Step 7
    --------|--------|--------|--------|--------|--------|--------|--------
          1 |      1 |      2 |      3 |      4 |      5 |      6 |
          2 |      7 |      8 |      9 |     10 |     11 |     12 |
          3 |     13 |     14 |     15 |     16 |     17 |     18 |
          4 |     19 |     20 |     21 |     22 |     23 |     24 |     25
          5 |     26 |     27 |     28 |     29 |     30 |     31 |
          6 |     32 |     33 |     34 |     35 |     36 |     37 |
          7 |     38 |     39 |     40 |     41 |     42 |     43 |
          8 |     44 |     45 |     46 |     47 |     48 |     49 |     50

To facilitate the transitions from one power to the next, which may require
additional move tuning (since moves should become bolder for powers closer to
0), a stone-specific pre-burnin (2) is performed. No samples are collected
during this stage, and another separate generation counter is used. By default,
this pre-burnin doubles the specified stone length, so if the user specifies
`.run(generations=1000)`, it will consist of 1000 generations for each stone.
This value can be changed by setting the `preburninGenerations` argument of the
`.run()` method.

Additionally, samples collected during the transition from one likelihood power
to another can be discarded as a stone-specific burnin (3), which is set to 25%
of the specified stone length by default. This value can be changed by setting
the `burninFraction` argument of the `.run()` method. The burnin phase shares
its generation counter with the sampling phase, so in the example above, the
first logged iteration would actually be 250 (with iterations 0--240 discarded
as burnin).

Finally, during the stone-specific sampling phase (4), samples are drawn from
a given power posterior. Its length is jointly determined by the `generations`
and `burninFraction` arguments of the `.run()` method, so in the example above,
it would consist of 750 generations for each stone (1000 specified generations
minus the 25% reserved for burnin).

In addition to the `.run()` method, which distributes stones among workers 
automatically, `powerPosterior` also supports more flexible parallelization
strategies via the `.runOneStone()` method, which allows the user to execute
a separate analysis for each power posterior. This increased flexibility comes
at the cost of lower-quality starting values: as `.runOneStone(index=<x>, ...)`
is not guaranteed to have been called before `.runOneStone(index=<x+1>, ...)`, 
the latter cannot be initialized with the output of the former, and will simply 
start from whatever state the sampler happens to currently hold.

Checkpointing

Stages 1 (global pre-burnin), 3 (stone-specific burnin), and 4 (stone-specific
sampling phase) can be checkpointed and resumed from a checkpoint. For highly
parallelized analyses where each stone is assigned to a separate worker, the
amount of time spent in stages 2--4 may be negligible compared to the amount
of time required to converge to the posterior. In such a case, it may be useful
to checkpoint and resume the global pre-burnin as follows:

    pow_p.burnin(generations=10000, tuningInterval=100,
                 checkpointFile="pow_p_burnin.ckp", checkpointInterval=100)
    pow_p.initializeFromCheckpoint("pow_p_burnin.ckp")
    pow_p.burnin(generations=10000, tuningInterval=100)

Note that the necessary checkpoint files can also be produced by regular MCMC
or MCMCMC analyses, so the `powerPosterior` sampler makes it possible to re-use
the output of such analyses without requiring a lengthy global pre-burnin stage
of its own.

Individual stones can also be checkpointed and resumed from a checkpoint. To do
so, the `.initializeFromCheckpoint()` method takes a `stones` argument that
allows the user to specify exactly which power posteriors should be restarted.
This can be a single index or a vector thereof. Following on the example above,
if the original analysis was interrupted during the very last step (step 7), we
may resume it as follows:

    pow_p.run(generations=1000, checkpointFile="analysis.ckp",
              checkpointInterval=100)
    pow_p.initializeFromCheckpoint("analysis.ckp", [25, 50])
    pow_p.run(generations=800)
    
When applied to individual stones, the checkpointing mechanism also records the
planned burnin length. Therefore, if the original analysis was interrupted
at generation 200 (i.e., still within the 250-generation-long burnin stage),
the commands above will cause the resumed analysis to complete the remaining 50
generations of burnin, followed by sampling for 750 generations (numbered 250
through 1000) -- exactly as though no interruption had occurred in the first
place. The `burninFraction` argument of the `.run()` call (which would normally
cause the first 0.25 times 800 = 200 generations to be discarded as burnin) is
ignored in this case, and a warning is printed to this effect.

If the analysis was interrupted at an earlier step, we may have to resurrect
entire sequences of stones. In this case, a nested vector (vector of vectors)
of stone indices is passed to `.initializeFromCheckpoint()` via its `stones`
argument. Each inner vector then corresponds to a sequence of stones to be
executed by one parallel worker. Only the first stone from each sequence needs
to have a checkpoint associated with it, as each subsequent stone will instead
use the output of the previous stone as its starting state. For example, if the
original analysis was interrupted during step 2, we could resume it as follows:

    pow_p.initializeFromCheckpoint("analysis.ckp",
              [ 2:6, 8:12, 14:18, 20:25, 27:31, 33:37, 39:43, 45:50 ])
    pow_p.run(generations=800)
    
Or as follows if it was interrupted during step 6:

    pow_p.initializeFromCheckpoint("analysis.ckp",
              [ [6], [12], [18], 24:25, [31], [37], [43], 49:50 ])

Note that in this case, the `.run()` call honors the assignment of stones to
workers specified in `.initializeFromCheckpoint()`, instead of automatically
generating an assignment of its own. An exception is thrown if the length of
the outer vector (i.e., the number of parallel workers requested by the user
for the resumed analysis) exceeds the number of workers available. Note also
that when the `stones` argument is specified using the range-based notation,
it makes a difference whether the range is enclosed in square brackets. The
following:

    pow_p.initializeFromCheckpoint( "analysis.ckp", 1:5 )
    
means "resurrect stones 1 through 5, and leave it up to the `.run()` method to
decide how to assign them to available workers". Checkpoint files must be
available for all 5 specified stones. In contrast,

    pow_p.initializeFromCheckpoint( "analysis.ckp", [1:5] )
    
means "resurrect stone 1, and then run stones 2 through 5 on the same parallel
worker, using the last sample from the previous stone as the starting point for
the next". In this case, we only require checkpoint files for stone 1.

Finally, checkpointing is also available for the `.runOneStone()` method:

    pow_p.runOneStone(index=16, generations=500, checkpointFile="analysis.ckp",
                      checkpointInterval=50)
    pow_p.initializeFromCheckpoint("analysis.ckp", 16)
    pow_p.runOneStone(index=16, generations=500)

## authors
Sebastian Hoehna
Michael Landis
John Huelsenbeck
David Černý
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
- citation: Lartillot N, Philippe H (2006). Computing Bayes factors using thermodynamic integration. Systematic Biology, 55(2):195-207.
  doi: 10.1080/10635150500433722
  url: https://academic.oup.com/sysbio/article-abstract/55/2/195/1620800
- citation: Xie W, Lewis PO, Fan Y, Kuo L, Chen M-H (2010). Improving marginal likelihood estimation for Bayesian phylogenetic model selection. Systematic Biology, 60(2):150-160.
  doi: 10.1093/sysbio/syq085
  url: https://academic.oup.com/sysbio/article-abstract/60/2/150/2461669
