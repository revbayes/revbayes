## name

mvRandomIntegerWalk

## title
Random walk on integers

## description
A move that performs random walk on an integer variable. The displacement of the random walk is exactly one step, either positive or negative.

## details

## authors
Sebastian Hoehna

## see_also
mvRandomNaturalWalk
mvRandomGeometricWalk

## example

p <- 0.8
x ~ dnGeom(p)

moves[1] = mvRandomIntegerWalk(x, weight=1.0)
monitors[1] = mvScreen(printgen=1000, x)

mymodel = model(p)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=20000,tuningInterval=100)
mymcmc.run(generations=200000)

## references
