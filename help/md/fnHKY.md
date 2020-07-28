## name
fnHKY

## title
The Hasegawa-Kishino-Yano (1985) nucleotide rate matrix

## description
The HKY85 model.

## details

## authors

## see_also

## example
        kappa ~ dnLognormal(0,1)
        pi ~ dnDirichlet( v(1,1,1,1) )
        Q := fnHKY(kappa,pi)

## references
- citation: Hasegawa, M. et al. Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. Journal of molecular evolution (1985) 22 (2): 160-174.
  doi: https://doi.org/10.1007/BF02101694
  url: https://link.springer.com/article/10.1007%2FBF02101694
