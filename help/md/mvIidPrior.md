## name
mvIidPrio
## title
Move to propose from prior
## description
This move proposes new values drawn from the prior.
## details
Using this move, one actually gets an independence sampler as the proposal doesn't depend on the current state. The move calls redraw based on the distribution attached to the random variable.
## authors
Sebastian Hoehna
## see_also
## example
x ~ dnUnif(0,10000)
moves[1] = mvIidPrior(x, weight=1.0)
monitors[1] = screenmonitor(printgen=1000, x)
mymodel = model(x)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.run(generations=200000)
## references
