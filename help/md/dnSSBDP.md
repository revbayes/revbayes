## name
dnSSBDP
## title
Sampled-speciation birth-death process
## description
Uses a data augmentation approach to sample hidden speciation events on
a birth-death tree.
## details
The sampled-speciation process is intended to be used in conjunction with
dispersal-extinction-cladogenesis (DEC) biogeographic models.
## authors
Michael Landis
## see_also
dnCDBDP
dnPhyloCTMCClado
## example
    # draw basic process parameters
    taxa <- [taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E")]
    root_age ~ dnUniform(0, 2)
    lambda ~ dnUniform(0, 1)
    mu ~ dnUniform(0, 1)
    rho <- 1

    # simulate tree
    tree ~ dnSSBDP(rootAge = root_age,
                   lambda  = lambda,
                   mu      = mu,           
                   rho     = rho,
                   taxa    = taxa)
## references
