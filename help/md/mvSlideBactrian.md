## name
mvSlideBactrian
## title
Bactrian Slide Move for MCMC 
## description
mvSlideBactrian is an MCMC move that proposes new values differently than the usual slide move
## details
Instead of small steps, it jumps further using a Bactrian distribution, which helps the chain explore the parameter space faster and reduces autocorrelation.
## authors
## see_also
mvSlide
mvScale
## example
x ~ dnUniform(0, 10)
moves.append(mvSlideBactrian(x, sigma=0.5, m=0.95, tune=TRUE, weight=1))
## references
