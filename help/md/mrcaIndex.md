## name
mrcaIndex
## title
## description
Returns the index of the most recent common ancestor (MRCA) of a given `clade`
in a given `tree`.
## details
## authors
Michael Landis
## see_also
tmrca
## example
    # Read in a tree
    tr <- readTrees(text="(((Lemur,Galago),Tarsius),((Pan,Homo),Alouatta));")[1]
    
    # Find the node index of the MRCA of Pan and Homo
    mrcaIndex(tr, clade("Homo", "Pan"))
    
## references
