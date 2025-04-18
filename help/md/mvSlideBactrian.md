## name
mvSlideBactrian
## title
Bactrian Slide Move for MCMC 
## description
mvSlideBactrian is an MCMC move that proposes new values differently than the usual slide move
## details
This move adds a random increment to the current value. The increment is drawn from a Bactrian distribution, which is symmetric around the current value but has low density at the center, creating a bimodal shape. As a result, it is less likely to propose very small changes, encouraging larger moves. This improves mixing efficiency and helps reduce autocorrelation in the MCMC chain
## authors
## see_also
mvSlide
mvScale
## example
x ~ dnUniform(0, 10)
moves.append(mvSlideBactrian(x, sigma=0.5, m=0.95, tune=TRUE, weight=1))
## references
- Citation: Z. Yang, & C.E. Rodríguez, Searching for efficient Markov chain Monte Carlo proposal kernels, Proc. Natl. Acad. Sci. U.S.A. 110 (48) 
  19307-19312, https://doi.org/10.1073/pnas.1311790110 (2013).