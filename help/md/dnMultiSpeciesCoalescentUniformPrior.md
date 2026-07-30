## name
dnMultiSpeciesCoalescentUniformPrior
## title
Multispecies coalescent distribution with uniform prior on effective population sizes
## description
Multispecies coalescent distribution describing how gene trees can be generated from within a species tree given effective population sizes. Requires an ultrametric species tree, an upper bound for the uniform prior on effective population sizes (a single real positive), and taxa with species and individual names.
## details
The species tree must be ultrametric.
This distribution uses a uniform prior on effective population sizes. As a consequence, effective population sizes are analytically integrated out and treated as nuisance parameters (Hey & Nielsen 2007).
## authors
Sebastian Hoehna
Bastien Boussau
## see_also
dnMultiSpeciesCoalescent
dnMultiSpeciesCoalescentInverseGamma
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
    
    # Let's pick a maximum effective population size of 50:
    popSize <- 50
    
    # Let's simulate gene trees now.
    # Taxon names:
    for (g in 1:n_genes) {
        for (i in 1:n_species) {
            for (j in 1:n_alleles) {
                taxa[g][(i-1)*n_alleles+j] <- taxon(taxonName="Species_"+i+"_"+j, speciesName="Species_"+i)
            }
        }
        geneTrees[g] ~ dnMultiSpeciesCoalescentUniformPrior(speciesTree=spTree, max=popSize, taxa=taxa[g])
        print(geneTrees[g])
    }
    
    # We can save the species tree:
    write(spTree, filename=dataFolder+"speciesTree")
    
    # Saving the gene trees:
    for (i in 1:(n_genes)) {
        write(geneTrees[i], filename=dataFolder+"geneTree_"+i+".tree")
    }

## references
- citation: Rannala B, Yang Z (2003). Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci. Genetics, 164(4):1645-1656.
  doi: 10.1093/genetics/164.4.1645
  url: https://academic.oup.com/genetics/article-abstract/164/4/1645/6050371
- citation: Heled J, Drummond AJ (2010). Bayesian inference of species trees from multilocus data. Molecular Biology and Evolution, 27(3):570-580.
  doi: 10.1093/molbev/msp274
  url: https://academic.oup.com/mbe/article/27/3/570/999753
- citation: Hey J, Nielsen R (2007). Integration within the Felsenstein equation for improved Markov chain Monte Carlo methods in population genetics. Proceedings of the National Academy of Sciences of the USA, 104(8):2785-2790.
  doi: 10.1073/pnas.0611164104
  url: https://www.pnas.org/content/104/8/2785
