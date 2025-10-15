## name
mvSpeciesNarrow
## title
Narrow-exchange joint move on species tree and gene trees for multispecies coalescent models.
## description
Makes a narrow-exchange move both in the species tree and in the gene trees that contain nodes of the relevant populations.
## details
The species tree must be ultrametric.

All the gene trees that evolved along the species tree according to some form of multispecies coalescent must be added to the move using the addGeneTreeVariable method.

This move jointly performs narrow exchange moves (Nearest-Neighbor Interchanges without branch length alterations) on the species tree and on gene trees, all of which must be ultrametric.

## authors
Sebastian Hoehna, Bastien Boussau
## see_also
mvSpeciesSubtreeScale
mvSpeciesSubtreeScaleBeta
mvSpeciesNodeTimeSlideUniform
mvSpeciesTreeScale
## example
    # We are going to save the trees we simulate in the folder simulatedTrees:
    dataFolder = "simulatedTrees/"
    
    # Let's simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:
    n_species <- 10
    n_genes <- 2
    n_alleles <- 3
    
    # We simulate an ultrametric species tree.
    # Species names:
    for (i in 1:n_species) {
        species[i] <- taxon(taxonName="Species_"+i, speciesName="Species_"+i)
    }
    spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy="uniform", condition="nTaxa", taxa=species)
    print(spTree)
    
    # Let's pick a constant effective population size of 50:
    popSize <- 50
    
    # Let's simulate gene trees now.
    # Taxon names:
    for (g in 1:n_genes) {
        for (i in 1:n_species) {
            for (j in 1:n_alleles) {
                taxa[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
            }
        }
        geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxa[g])
        print(geneTrees[g])
    }
    
    # Set my move index:
    mi = 0
    move_species_narrow_exchange = mvSpeciesNarrow( speciesTree=spTree, weight=5 )
    for (i in 1:n_genes) {
        move_species_narrow_exchange.addGeneTreeVariable( geneTrees[i] )
    }
    moves[++mi] = move_species_narrow_exchange
    
    # We get a handle on our model.
    # We can use any node of our model as a handle; here we choose to use the topology.
    mymodel = model(spTree)
    
    # Monitors to check the progression of the program:
    monitors[1] = mnScreen(printgen=10, spTree)
    
    # Here we use a plain MCMC. You could also use mcmcmc with heated chains.
    mymcmc = mcmc(mymodel, monitors, moves, nruns=4)
    mymcmc.run(generations=1000)
    mymcmc.operatorSummary()
    
## references
- citation: Hoehna S, Drummond AJ (2012). Guided tree topology proposals for Bayesian phylogenetic inference. Systematic Biology 61(1):1-11.
  doi: 10.1093/sysbio/syr074
  url: https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syr074
- citation: Jones G (2016). Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Journal of Mathematical Biology, 74:447-467.
  doi: 10.1007/s00285-016-1034-0
  url: http://www.indriid.com/2016/2016-06-01-STACEY.pdf
