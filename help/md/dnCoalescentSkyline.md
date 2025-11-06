## name
dnCoalescentSkyline
## title
Heterochronous and homochronous skyline coalescent process
## description
The skyline coalescent process specifies a probability density on genealogies, both node ages and the topology. It is used for both heterochronous samples and homochronous samples.
## details
The underlying theory of the skyline coalescent implemented here is Kingman's coalescent. The implementation here assumes haploid individuals, so for diploid study systems one needs to multiply the effective population size by 2 and the true effective population size in units of individuals needs to be divided by 2 afterwards.

The coalescent process is parameterized with `theta`, which here stands for the effective population size (not 4 * Ne * mu). For detailed examples see https://revbayes.github.io/tutorials/coalescent/
## authors
Ronja Billenstein
Andrew Magee
Sebastian Hoehna
## see_also
dnCoalescent
dnCoalescentDemography
## example
    NUM_INTERVALS = ceil(n_taxa / 5)
    for (i in 1:NUM_INTERVALS) {
        
        pop_size[i] ~ dnUniform(0,1E6)
        pop_size[i].setValue(100.0)
        moves.append( mvScale(pop_size[i], lambda=0.1, tune=true, weight=2.0) )
    
    }
    
    # next we specify a prior on the number of events per interval
    # we use a multinomial prior offset to have at least one event per interval
    # first, specify the offset
    num_events_pi <- rep(1, NUM_INTERVALS)
        
    # next, specify the prior for the multinomial distribution
    num_e_simplex_init <- rep(1, NUM_INTERVALS)
    num_e_simplex <- simplex(num_e_simplex_init)
    
    # calculate the number of coalescent events that we distribute over the intervals
    n_multi <- n_taxa-1-NUM_INTERVALS
    
    # draw the coalescent events into intervals
    number_events_pi ~ dnMultinomial(p=num_e_simplex, size=n_multi)
    
    # compute the actual number of events per interval, so the drawn number plus offset
    final_number_events_pi := num_events_pi + number_events_pi
    
    moves.append( mvIidPrior(x=number_events_pi) )
    
    # the time tree is a stochastic node modeled by the constant-rate coalescent process (dnCoalescent)
    psi ~ dnCoalescentSkyline(theta=pop_size, events_per_interval=final_number_events_pi, method="events", taxa=taxa)
    
    interval_times := psi.getIntervalAges()
    
    root_height := psi.rootAge()

    # continue as usual to either clamp the genealogy or infer the genealogy based on sequence data

## references
- citation: Billenstein R, Höhna S (2024). Comparison of Bayesian coalescent skyline plot models for inferring demographic histories. Molecular Biology and Evolution, 41(5):msae073.
  doi: 10.1093/molbev/msae073
  url: https://academic.oup.com/mbe/article/41/5/msae073/7648822
