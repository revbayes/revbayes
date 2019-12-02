## name
mvSpeciesSubtreeScale
## title
Subtree scale move on species tree and gene trees for multispecies coalescent models.
## description
Makes a subtree scale move both in the species tree and in the gene trees that contain nodes of the relevant populations. Tree topologies are not altered.
## details
The species tree must be ultrametric.

All the gene trees that evolved along the species tree according to some form of multispecies coalescent must be added to the move using the addGeneTreeVariable method.

This move jointly performs a subtree scale move (a whole subtree is scaled up or down, keeping the topology fixed) on the species tree and on gene trees, all of which must be ultrametric.

How this works: we pick a random node which is not the root.
Then, we uniformly pick an age between the parent and the oldest sampled descendant.
The picked subtree is then scaled to this new age.
All gene-trees that are present in the population will be scaled accordingly.

## authors
Sebastian Hoehna, Bastien Boussau
## see_also
mvSpeciesNodeTimeSlideUniform
mvSpeciesSubtreeScaleBeta
mvSpeciesNarrow
mvSpeciesTreeScale
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
	# let's pick a constant effective population size of 50:
	popSize <- 50
	# let's simulate gene trees now:
	# taxa names:
	for (g in 1:n_genes) {
	  for (i in 1:n_species) {
	    for (j in 1:n_alleles) {
	        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
	    }
	  }
	  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
	  print(geneTrees[g])
	}
	# We can save the species tree and the gene trees: 
	write(spTree, filename=dataFolder+"speciesTree")
	# Saving the gene trees
	for (i in 1:(n_genes)) {
	  write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
	}
	# set my move index
	mi = 0
	move_species_subtree_scale = mvSpeciesSubtreeScale( speciesTree=spTree, weight=5 )
	for (i in 1:n_genes) {
	   move_species_subtree_scale.addGeneTreeVariable( geneTrees[i] )
	}
	moves[++mi] = move_species_subtree_scale
	# We get a handle on our model.
	# We can use any node of our model as a handle, here we choose to use the topology.
	mymodel = model(spTree)
	# Monitors to check the progression of the program
	monitors[1] = mnScreen(printgen=10, spTree)
	# Here we use a plain MCMC. You could also set nruns=2 for a replicated analysis
	# or use mcmcmc with heated chains.
	mymcmc = mcmc(mymodel, monitors, moves, nruns=4)
	mymcmc.run(generations=1000)
	mymcmc.operatorSummary()
	
## references
- citation: "Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Sebastian\
    \ H\xF6hna, Alexei J. Drummond. Syst Biol (2012) 61 (1): 1-11."
  doi: https://doi.org/10.1093/sysbio/syr074
  url: https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syr074
- citation: Algorithmic improvements to species delimitation and phylogeny estimation
    under the multispecies coalescent. Graham Jones.  Journal of Mathematical Biology,
    2016.
  doi: https://doi.org/10.1007/s00285-016-1034-0
  url: http://link.springer.com/article/10.1007/s00285-016-1034-0
