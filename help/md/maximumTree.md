## name
maximumTree
## title
Maximum tree function to build a species tree.
## description
Builds the maximum species tree given several ultrametric gene trees.
## details
The maximum species tree is a consistent estimate of the species tree under the multispecies coalescent model, if the gene trees are correct and the effective population size constant along the species tree.
## authors
Bastien Boussau
## see_also
## example
	# We simulate a species tree and gene trees and reconstruct a species tree using maximum tree:
	# Let’s simulate a species tree with 10 taxa, 5 gene trees, 1 alleles per species:
	n_species <- 10
	n_genes <- 5
	n_alleles <- 2
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
	        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
	    }
	  }
	  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])
	}
	# Let's compute the maximum tree: 
	recTree <- maximumTree(geneTrees)
	print(recTree)
	
## references
- citation: High-resolution species trees without concatenation. Scott V. Edwards,
    Liang Liu, and Dennis K. Pearl . PNAS April 3, 2007 vol. 104 no. 14 .
  doi: null
  url: http://www.pnas.org/content/104/14/5936.full
- citation: 'Maximum tree: a consistent estimator of the species tree. Liu L, Yu L,
    Pearl DK.  Journal of Mathematical Biology, 2010. Jan;60(1):95-106.'
  doi: https://doi.org/10.1007/s00285-009-0260-0
  url: https://link.springer.com/article/10.1007%2Fs00285-009-0260-0
