## name
mvScaleBactrian
## title
Scaling Move Employing Bactrian Distribution
## description
Scales a parameter by a Bactrian-distributed factor. 
## details
Proposes multiplicative changes to a real-valued parameter using a Bactrian
kernel -- a bimodal distribution centered at zero, obtained as a mixture of two
unimodal component distributions. Specifically, `mvScaleBactrian` scales
parameter values by a random factor of exp(lambda * delta), where lambda is a
tuning parameter and delta is drawn from a mixture of two normal distributions
following Yang & Rodríguez (2013: Supplementary Information, Eq. 19): 

    u ~ Unif(0, 1)
    x ~ N(0, 1)
    delta = m + x * sqrt(1 - m^2)  if u < 0.5
          = -m + x * sqrt(1 - m^2) otherwise

with m set to 0.95. As a result, the move is less likely to propose very small
steps and encourages larger changes, which improves mixing efficiency and
reduces autocorrelation. Since multiplicative updates do not change the sign,
`mvBactrianScale` is most often applied to parameters that are constrained to
be positive.

## authors
## see_also
mvScale
mvScaleBactrianCauchy
## example
    moves = VectorMoves()
    speciation_rate ~ dnGamma(2, 4)
    moves.append( mvScaleBactrian(speciation_rate, weight=5) )

## references
- citation: Yang Z, Rodríguez CE (2013). Searching for efficient Markov chain Monte Carlo proposal kernels. Proc. Natl. Acad. Sci. USA, 110(48):19307-19312.
  doi: 10.1073/pnas.1311790110
  url: https://www.pnas.org/doi/full/10.1073/pnas.1311790110
