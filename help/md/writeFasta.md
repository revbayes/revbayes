## name
writeFasta
## title
FASTA file writing function
## description
This function writes out a FASTA formatted file given 
data of class `AbstractHomologousDiscreteCharacterData`.
Filename is specified using the `fn` argument.
## details
## authors
## see_also
`writeNexus`, `writeCharacterDataDelimited`
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
    # this will print a FASTA file with a simulated molecular matrix
    # to the working directory
    writeFasta(filename="test.fasta", full_model)
    
## references
- citation: Pearson, William R., and David J. Lipman. "Improved tools for biological sequence comparison." Proceedings of the National Academy of Sciences 85.8 (1988): 2444-2448.
- url: https://www.pnas.org/content/85/8/2444.short
- doi: https://doi.org/10.1073/pnas.85.8.2444
