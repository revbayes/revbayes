## name
Set value
## title
Set the value of a stochastic variable

## description
`x.setValue(value)` sets the value of the stochastic variable `x` to `value`.

## details
`x.setValue()` allows calculations to be evaluated at a given value of `x`, whilst allowing the value
of `x` to be redrawn, or to vary during MCMC.

`.setValue()` allows an MCMC run to be initialized with plausible values,
which can expedite convergence when priors are broad
([example](https://revbayes.github.io/tutorials/divrate/branch_specific.html#specifying-the-model)),
and can be useful when [debugging MCMC runs](https://revbayes.github.io/tutorials/mcmc_troubleshooting/#starting-values).

## authors

## see_also
clamp
## example
x ~ dnNormal(1, 1)

# Evaluate P(x) at x = 1
x.setValue(1)
x.probability()

# Modify the observed value of x
x.redraw()
x.probability()

# Initialize an MCMC run with a specific value
x.setValue(40000)
mcmc(model(x), [mnScreen(x)], [mvScale(x)]).run(generations = 5)


## references
