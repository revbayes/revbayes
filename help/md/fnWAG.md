## name
fnWAG
## title
WAG (Whelan and Goldman) Empirical Amino Acid Substitution Rate Matrix 
## description
The WAG model generates a rate matrix based on the Whelan and Goldman (WAG) substitution model for amino acid evolution

## details
The WAG model, a empirical model of amino acid replacement derived using an approximate maximum-likelihood method from 3,905 sequences across 182 protein families. It outperforms widely used models like Dayhoff and JTT in terms of accuracy and likelihood for phylogenetic analysis, aiming to provide better evolutionary tree estimates and applications in sequence alignment, database searches, and protein structure prediction.
## authors

## see_also
fnDayhoff
fnJones
fnLG
## example
  #WAG model with estimated frequencies 
  pi ~ dnDirichlet( rep(1,20) )
  Q := fnWAG(pi)

  #WAG model with fixed frequencies
  Q2 := fnWAG()



## references
- citation: Simon Whelan, Nick Goldman, A General Empirical Model of Protein Evolution Derived from Multiple Protein Families Using a Maximum-Likelihood Approach, Molecular Biology and Evolution, Volume 18, Issue 5, May 2001, Pages 691–699,
  url: https://doi.org/10.1093/oxfordjournals.molbev.a003851 
