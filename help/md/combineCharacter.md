## name
combineCharacter
## title
Concatenate character matrices
## description
Creates a new data matrix by concatentating the provided data matrices (by order).
## details
## authors
Sebastian Hoehna
## see_also
## example
	# read in character data for locus_1
	locus_1 = readContinuousCharacterData("locus_1.nex")
	# read in character data for locus_2
	locus_2 = readContinuousCharacterData("locus_2.nex")
	# create concated locus for 1+2 (union of taxa)
	locus_1_and_2 = concatenate( locus_1, locus_2 )
	
## references
