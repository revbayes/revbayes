## name
mvLayeredScaleProposal
## title
Rescales all the subtrees below some age.
## description
Makes a subtree scale move on all subtrees below a given age in the tree. Tree topology is not altered.
## details
The tree must be ultrametric.

An age is randomly drawn between the root age and the age of the oldest tip. Then all subtrees below this age are scaled up or down depending on a scaler drawn from an exponential distribution.

## authors
Bastien Boussau
## see_also
mvSubTreeScale
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
    move_species_subtree_scale = mvLayeredScaleProposal( speciesTree=spTree, weight=5 )
    for (i in 1:n_genes) {
        move_species_subtree_scale.addGeneTreeVariable( geneTrees[i] )
    }
    moves[++mi] = move_species_subtree_scale
    
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
