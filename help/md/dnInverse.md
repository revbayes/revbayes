## name
dnInverse

## title
Inverse distribution 

## description
`dnInverse()` inverts a probability distribution.

`dnInverse(x).probability()` returns `1 / x.probability()`; 
`dnInverse(x).lnProbability()` returns `-x.lnProbability()`.

This provides a way to perform inference using conditional probabilities,
for example where 
Pr(x | Model, Condition) = Pr(x | Model) / Pr(Condition is satisfied)

## details
## authors
Martin R. Smith

## see_also
## example
```
# Compute Pr(x = 1 | x ~ exp(y), y = 1)

# First we define the distributions from which x and y are drawn

# y may take the values 1 or 2 with probabilities 0.4, 0.6
p_y := simplex(0.4, 0.6)
y ~ dnCategorical( p_y )
inv_y ~ dnInverse(dnCategorical( p_y ))

y.clamp(1)
inv_y.clamp(1)

y.probability()      # 0.4 = 2 / 5
inv_y.probability()  # 2.5 = 5 / 2

# Compute the joint probability Pr(x, y)
function PrXandY (x_value, y_value) {
    x_given_y ~ dnExponential( y_value )
    x_given_y.clamp(x_value) # To calculate Pr(x | y)

    y.clamp(y_value) # To calculate Pr(y)

    return(x_given_y.probability() * y.probability())
}

PrXandY(1, 1) # Pr(x = 1, y = 1)

# If we wish to calculate likelihood conditioned on y = 1, we need to compute
# Pr(x = 1 | y = 1)

# Here it is trivial to compute this directly, as in PrXandY, but in more complex
# cases it may be easier to use Pr(x = 1 | y = 1) := Pr(x = 1, y = 1) / Pr(y = 1)

PrXandY(1, 1) * inv_y.probability()
```

## references
