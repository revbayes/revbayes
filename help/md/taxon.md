## name
taxon
## title
Taxon object
## description
The taxon function creates a Taxon object.
## details
Each Taxon object records that taxon's name in addition to other information, such as its age range.
A taxon constructed without ages has no age range. Use `hasAgeRange()` to distinguish this from
an explicitly supplied range, including [0,0]. Character matrices and trees contain
Taxon vectors (Taxon[]) that are used to match leaf nodes to data entries for phylogenetic analyses.
For multispecies coalescent analyses, Taxon objects are also used to assign species memberships to
individuals.

Use `setAgeRange(minAge, maxAge)` to replace both age bounds together. The bounds must be
nonnegative, must not be NaN, and must satisfy `minAge <= maxAge`; equal bounds specify an exact age.
An invalid call leaves the existing range unchanged.

`getMinAge()` and `getMaxAge()` return the endpoints of an existing range; `getAge()` is an alias
for `getMinAge()`. These getters throw if no range has been supplied. Use `hasAgeRange()` to check
before querying a potentially absent range.
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

	# Set both bounds of a fossil's age range together
	fossil = taxon("Fossil_species")
	fossil.setAgeRange(minAge=5, maxAge=6)
	
## references
