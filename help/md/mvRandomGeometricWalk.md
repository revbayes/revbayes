## name
mvRandomGeometricWalk

## title
Geometric random walk

## description
A move that performs geometric random walk on an integer variable. The displacement of the random walk is drawn from a geometric distribution, mirrored for positive and negative steps.

## details

## authors
Sebastian Hoehna

## see_also
mvRandomNaturalWalk
mvRandomIntegerWalk

## example

p <- 0.8
x ~ dnGeom(p)

moves[1] = mvRandomGeometricWalk(x, weight=1.0)
monitors[1] = mvScreen(printgen=1000, x)

mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000)

## references
