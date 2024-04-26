## name
dnMultiSpeciesCoalescentInverseGamma
## title
Multispecies coalescent Distribution with inverse gamma prior on effective population sizes
## description
Multispecies coalescent distribution describing how gene trees can be generated from within a species tree given effective population sizes. Requires an ultrametric species tree, parameters of an inverse gamma prior on effective population sizes, and taxa with species and individual names.
## details
The species tree must be ultrametric.
Parameters of an inverse gamma prior on effective population sizes must be provided.
This distribution uses a conjugate prior on effective population sizes. As a consequence, effective population sizes are integrated out and treated as nuisance parameters.

If you are interested in reconstructing ancestral effective population sizes, use dnMultiSpeciesCoalescent.

## authors
Sebastian Hoehna, Bastien Boussau
## see_also
dnMultiSpeciesCoalescent
dnMultiSpeciesCoalescentUniformPrior
## example
	# We are going to save the trees we simulate in the folder simulatedTrees:
	dataFolder = "simulatedTrees/" 
	# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
	n_species <- 10
	n_genes <- 2
	n_alleles <- 3
	# we simulate an ultrametric species tree:
	# Species names:
	for (i in 1:n_species) {
	        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
	}
	spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
	print(spTree)
	# let's pick constant parameters for the inverse gamma distribution:
	alpha <- 3
	beta <- 0.003
	# let's simulate gene trees now:
	# taxa names:
	for (g in 1:n_genes) {
	  for (i in 1:n_species) {
	    for (j in 1:n_alleles) {
	        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
	    }
	  }
	  geneTrees[g] ~ dnMultiSpeciesCoalescentInverseGamma(speciesTree=spTree, shape=alpha, scale=beta, taxa=taxons[g])
	  print(geneTrees[g])
	}
	# We can save the species tree and the gene trees: 
	write(spTree, filename=dataFolder+"speciesTree")
	# Saving the gene trees
	for (i in 1:(n_genes)) {
	  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
	}
	
## references
- citation: ' Algorithmic improvements to species delimitation and phylogeny estimation
    under the multispecies coalescent. Jones G.  Journal of Mathematical Biology.
    2016.'
  doi: 'DOI: 10.1007/s00285-016-1034-0'
  url: http://www.indriid.com/2016/2016-06-01-STACEY.pdf
