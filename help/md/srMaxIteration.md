## name
srMaxIteration
## title
Maximum iteration stopping rule
## description
Terminates an MCMC run once the specified number of iterations have been
performed.
## details
For an analysis conducted from scratch, employing the maximum iteration rule is
exactly equivalent to specifying the `generations` argument:

    mcmc.run( [srMaxIteration(10000)] )

as well as the following:

    mcmc.run(generations=10000)

will both stop the analysis once it has completed 10,000 generations. However,
the two `.run()` calls above will behave differently when they follow a call to
`.initializeFromCheckpoint()` -- i.e., when resuming a checkpointed analysis.
For example, if such an analysis was interrupted at generation 5000, resuming
it with `generations=10000` will cause it to be extended by another 10,000
generations, yielding a total length of 15,000 generations. In contrast,
specifying `srMaxIteration(10000)` upon resumption will stop the analysis once
it has reached a total length of 10,000 generations, limiting the resumed run
to just 5000 generations.
## authors
## see_also
mcmc
srMaxTime
## example
    # Binomial example: estimate success probability given 7 successes out of 20 trials
    r ~ dnExp(10)
    p := Probability(ifelse(r < 1, r, 1))
    n <- 20
    k ~ dnBinomial(n, p)
    k.clamp(7)
    mymodel = model(k)
    
    moves = VectorMoves()
    moves.append( mvSlide(r, delta=0.1, weight=1) )
    
    paramFile = "parameters.log"
    
    monitors = VectorMonitors()
    monitors.append( mnModel(filename=paramFile, printgen=100, p) )
    
    # Stop when 1000 iterations have been completed
    stopping_rules[1] = srMaxIteration(1000)
    
    # Create the MCMC object
    mymcmc = mcmc(mymodel, monitors, moves)
    
    # Begin the MCMC run
    mymcmc.run(rules = stopping_rules)

## references
