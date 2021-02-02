## name
clade
## title
Clade
## description
A clade is a subtree within a phylogeny.
## details
Clades are defined in terms of a taxon set and a shared tree topology. In phylogenetic analyses, clades are generally used (a) to constrain tree topologies to match provided taxon relationships, (b) to identify the most recent common ancestor of a taxon set within a phylogeny, or (c) to apply node age calibrations on particular nodes in the phylogeny.
## authors
Will Pett
Sebastian Hoehna
Michael Landis
## see_also
dnConstrainedTopology
tmrca
mrcaIndex
## example
	# read in a tree
	phy = readTrees("primates.tre")[1]
	# get taxa from the tree
	taxa = phy.taxa()
	# create a clade for (1,2) using taxon objects
	clade_12 = clade( taxa[1], taxa[2] )
	# create a clade for (1,2,3)
	clade_123 = clade( taxa[3], clade_12 )
	# create a clade for (4,5) using taxon names
	clade_45 = clade( "Taxon_4", "Taxon_5" )
	# create a negative clade constraint
	clade_not_23 = clade( taxa[2], taxa[3], negative=true )
	# create an optional clade constraint
	clade_67 = clade( taxa[6], taxa[7] )
	clade_68 = clade( taxa[6], taxa[8] )
	clade_67_or_68 = clade( clade_67, clade_68, optional_match=true )
	
## references
