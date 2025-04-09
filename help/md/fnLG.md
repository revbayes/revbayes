## name
fnLG
## title
LG (Le and Gascuel) Amino Acid Substitution Rate Matrix 
## description
Generates a rate matrix based on the Le and Gascuel (LG) substitution model for amino acid evolution. 
## details
The LG model is an empirical model of amino acid replacement derived from large-scale protein alignments. It refines substitution rate estimation by explicitly considering site-specific rate variation.
## authors
## see_also
fnDayhoff
fnJones
fnWAG
## example
  # LG model with estimated frequencies
  pi ~ dnDirichlet( rep(1,20) )
  Q := fnLG(pi)

  # LG model with fixed frequencies
  Q2 <- fnLG()

## references
- citation: Le SQ, Gascuel O (2008). An improved general amino acid replacement matrix. Molecular Biology and Evolution, 25(7):1307-1320. 
  doi: 10.1093/molbev/msn067
  url: https://academic.oup.com/mbe/article/25/7/1307/1041491
