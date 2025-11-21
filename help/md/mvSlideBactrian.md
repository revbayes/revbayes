## name
mvSlideBactrian
## title
Sliding-Window Mode Employing Bactrian Distribution
## description
Updates a parameter by a Bactrian-distributed increment.
## details
Proposes additive changes to a real-valued parameter using a Bactrian kernel
-- a bimodal distribution centered at zero, obtained as a mixture of two
unimodal component distributions. Specifically, `mvSlideBactrian` updates
the current value by adding a random increment of (lambda * delta), where
lambda is a tuning parameter and delta is drawn from a mixture of two normal
distributions following Yang & Rodríguez (2013: Supplementary Information,
Eq. 19): 

    u ~ Unif(0, 1)
    x ~ N(0, 1)
    delta = m + x * sqrt(1 - m^2)  if u < 0.5
          = -m + x * sqrt(1 - m^2) otherwise

with m set to 0.95. As a result, the move is less likely to propose very small
steps and encourages larger changes, which improves mixing efficiency and
reduces autocorrelation.

## authors
## see_also
mvSlide
## example
    moves = VectorMoves()
    x ~ dnNormal(0, 2)
    moves.append( mvSlideBactrian(x, tune=TRUE, weight=1) )

## references
- citation: Yang Z, Rodríguez CE (2013). Searching for efficient Markov chain Monte Carlo proposal kernels. Proc. Natl. Acad. Sci. USA, 110(48):19307-19312.
  doi: 10.1073/pnas.1311790110
  url: https://www.pnas.org/doi/full/10.1073/pnas.1311790110
