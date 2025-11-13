## name
writePhylip
## title
Phylip file writing function
## description
This function writes out a phylip formatted file given
data of class `AbstractHomologousDiscreteCharacterData`.
Filename is specified using the `filename` argument.
## details
## authors
## see_also
writeCharacterDataDelimited
writeFasta
writeNexus
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
    
    # this will print a Phylip file with a simulated molecular matrix
    # to the working directory
    writePhylip(filename="test.phy", full_model)
