## name
writeNexus
## title
Nexus file writer
## description
Function for writing a nexus file.  
## details
The first argument is the filename  to write to and this must be a string.
The second argument is a data object that must be some character matrix. 
This data matrix could be a morphological matrix, a molecular matrix, or a tree.

## authors
## see_also
`writeFasta`, `writeCharacterDataDelimited`, `write`
## example
    # let's make up some taxa
    taxa = v("horse", "whale", "unicorn", "narwhal")
    # simulate a tree
    tau ~ dnUniformTimeTree(rootAge=1, taxa=taxa)
    # we also need a molecular substitution model
    molecular_model := fnJC(4)
    # together these form a continuous time Markov chain over the tree
    full_model ~ dnPhyloCTMC(tree=tau, Q=molecular_model, nSites = 100, type="DNA")
    # this will print a Nexus file with a simulated molecular matrix
    # to the working directory
    writeNexus(filename="test.nex", full_model)
    
    
## references
- citation: David R. Maddison, David L. Swofford, Wayne P. Maddison, Nexus: An Extensible File Format for Systematic Information, Systematic Biology, Volume 46, Issue 4, December 1997, Pages 590–621,
  doi: https://doi.org/10.1093/sysbio/46.4.590
  url: https://academic.oup.com/sysbio/article/46/4/590/1629695