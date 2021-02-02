## name
tmrca
## title
Find the time to the most recent common ancestor
## description
Finds the most recent common ancestor (TMRCA) of a clade of taxa on a tree.
## details
## authors
## see_also
`clade`
## example
    # let's make up some taxa 
    taxa = v("horse", "whale", "unicorn", "narwhal")
    # convert these to the taxon datatype
    for(i in 1:4) { taxa[i] = taxon(taxa[i]) }
    # simulate a tree
    tau ~ dnUniformTimeTree(rootAge=1, taxa=taxa)
    # we also need a molecular substitution model
    molecular_model := fnJC(4)
    # together these form a continuous time Markov chain over the tree
    full_model ~ dnPhyloCTMC(tree=tau, Q=molecular_model, nSites = 100, type="DNA")
    # need to make a clade
    horned_animals <- clade(taxa[3], taxa[4])
    tmrca(tau, horned_animals)
## references
