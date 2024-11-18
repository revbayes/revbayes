## name
Clamp
## title
Clamp a stochastic variable to a fixed/observed value
## description
`x.clamp(data)` fixes the value of the stochastic variable `x` to the observation `data`, and marks the variable as corresponding to an observation.

## details
Once clamped, the value of `x` is thus expected to remain constant, unless `x` is subsequently unclamped – either explicitly with `x.unclamp()`, or implicitly with `x.clamp(different_data)`.

`x.setValue()` evaluates probabilities at a specific value of `x` without fixing the value.

## authors

## see_also
setValue
unclamp
## example
x ~ dnNormal(1, 1)
y ~ dnNormal(2, 2)
z := x * y

# Set the observed value of x
x.clamp(1)

# Evaluate P(z) at x = X, y = 1
y.setValue(1)
z.probability()

# Modify the observed value of x
x.clamp(2) # equivalent to x.unclamp(); x.clamp(2)
z.probability()

# Evaluate P(z) at a different value of y
y.redraw()
print(y)
z.probability()

# Because x is clamped, it is invalid to call x.redraw() or mvSlide(y)
# x will remain constant during MCMC, whereas y will be inferred.
mcmc(model(z), [mnScreen(x, y)], [mvSlide(y)]).run(generations = 5)


## references
