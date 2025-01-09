## name
dnCoalescent
## title
Constant population size Coalescent process
## description
The constant population size coalescent process specifies a probability density on genealogies, both node ages and the topology.
## details
The underlying theory of the constant population size Coalescent implemented here is Kingman's Coalescent. The implementation here assumes haploid individuals, so for diploid study systems one needs to multiply the effective population size by 2 and the true effective population size in units of individuals needs to be divided by 2 afterwards.
The Coalescent process is parameterized with theta, which here stands for the effective population size (not 4 * Ne * mu). For detailed examples see https://revbayes.github.io/tutorials/coalescent/
## authors
Ronja Billenstein
Andrew Magee
Sebastian Hoehna
## see_also
dnCoalescentSkyline
dnCoalescentDemography
## example

# specify a prior distribution on the constant population size
pop_size ~ dnUniform(0,1E6)
moves.append( mvScale(pop_size, lambda=0.1, tune=true, weight=2.0) )

# specify the coalescent process.
# note that you need to have a vector of taxa
psi ~ dnCoalescent(theta=pop_size, taxa=taxa)

# for monitoring purposes, you may want the root age
root_height := psi.rootAge()

# continue as usual to either clamp the genealogy or infer the genealogy based on sequence data

## references
- citation: Comparison of Bayesian Coalescent Skyline Plot Models for Inferring Demographic Histories. Billenstein, Ronja and Höhna, Sebastian (2024) Molecular Biology and Evolution, 41(5):msae073.
  doi: https://doi.org/10.1093/molbev/msae073
  url: https://academic.oup.com/mbe/article/41/5/msae073/7648822
