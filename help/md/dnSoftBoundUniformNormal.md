## name
dnSoftBoundUniformNormal
## title
Soft-bounded uniform-normal distribution
## description
A soft-bounded uniform distribution with one or two normally distributed tails
outside the domain of the uniform component.
## details
Draws from a uniform distribution defined between the `min` and `max` values
with probability `p`. The tail probability (1 - `p`) is assigned to values
smaller than `min` (if `boundary="lower"`), values larger than `max` (if
`boundary="upper"`), or apportioned equally between the lower and upper tails
(by default, when `boundary="both"`). The probability density outside the
[`min`, `max`] interval is that of a normal distribution whose standard
deviation can be either automatically computed from `p` or directly specified
using the `sd` argument. Exactly one of `sd` and `p` must be specified, with
the unspecified value automatically computed so that the density is continuous
at `min` and `max`.
## authors
Sebastian Hoehna
## see_also
dnNormal
dnUniform
## example
	# Create a simple model (unclamped)
    calib ~ dnSoftBoundUniformNormal( 13.6, 25, p=0.95, "upper" )
    mymodel = model(calib)
    
    # Create a move vector and a monitor vector
	moves[1] = mvSlide(calib, delta=0.1, weight=1.0)
	monitors[1] = mnScreen(calib, printgen=1000)
	
    # Use MCMC to draw samples from the specified distribution   
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.run(generations=200000)

## references
