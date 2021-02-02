## name
readTrees
## title
Function to read in trees.
## description
Reads trees from a file containing trees (Nexus, Phylip or Newick accepted), or from a string containing Newick representations of trees.
## details
Either a file name (with the file argument) or a string (with the text argument) must be provided as argument. If both are provided, trees will be read from both sources.
## authors
Bastien Boussau
## see_also
write
readBranchLengthTrees
readCharacterData
readCharacterDataDelimited
readContinuousCharacterData
readDiscreteCharacterData
readDataDelimitedFile
readCharacterData
## example
	trees=readTrees(text="(a,(b,c));\n(d:0.1,(e:0.1,f:0.1):0.1);")
	print(trees)
	
## references
