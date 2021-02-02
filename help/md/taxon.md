## name
taxon
## title
Taxon object
## description
The taxon function creates a Taxon object.
## details
Each Taxon object records that taxon's name in addition to other information, such as age (which is non-zero for fossils). Character matrices and trees contain Taxon vectors (Taxon[]) that are used to match leaf nodes to data entries for phylogenetic analyses. For multispecies coalescent analyses, Taxon objects are also used to assign species memberships to individuals.
## authors
Michael Landis
## see_also
readTaxonData
## example
	# we can create a Taxon object
	taxon_gorilla = taxon("Gorilla_gorilla")
	# we can create a dummy vector of Taxon objects for simulation
	for (i in 1:10) { taxa[i] = taxon("Taxon"+i) }
	phy ~ dnBDP(lambda=1, mu=0, rootAge=1, taxa=taxa)
	# retrieve the taxon list for 'phy'
	phy.taxa()
	
## references
