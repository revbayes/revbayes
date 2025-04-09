## name
fnWAG
## title
WAG (Whelan and Goldman) Amino Acid Substitution Rate Matrix 
## description
Generates a rate matrix based on the Whelan and Goldman (WAG) substitution model for amino acid evolution.
## details
The WAG model is an empirical model of amino acid replacement derived using an approximate maximum-likelihood method from 3,905 sequences across 182 protein families. It outperforms previous models like Dayhoff and JTT in terms of accuracy and likelihood for phylogenetic analysis, aiming to provide better evolutionary tree estimates and applications in sequence alignment, database searches, and protein structure prediction.
## authors
## see_also
fnDayhoff
fnJones
fnLG
## example
  # WAG model with estimated frequencies 
  pi ~ dnDirichlet( rep(1,20) )
  Q := fnWAG(pi)

  # WAG model with fixed frequencies
  Q2 <- fnWAG()

## references
- citation: Whelan S, Goldman N (2001). A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. Molecular Biology and Evolution, 18(5):691-699.
  doi: 10.1093/oxfordjournals.molbev.a003851
  url: https://academic.oup.com/mbe/article/18/5/691/1018653 
