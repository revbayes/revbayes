## name
dnPhyloBrownianProcessStateDependent (dnPhyloBMSD)
## title
Phylogenetic state-dependent Brownian process
## description
Univariate Brownian process over an augmented phylogeny (i.e., discrete character history)
## details
The phylogenetic state-dependent Brownian process is collapsed from the phylogenetic state-dependent Ornstein-Uhlenbeck process, in which `alpha` is always 0. The probability is computed using a pruing algorithm. Specifically, `dnPhyloBMSD` uses an augmented tree structure, i.e., the character history with different character states (also referred to as "regimes"; Hansen 1997). Shifts in the character state is discrete and does not occur concurrently with speciation events, and can be represented by nodes of degree 2. As such, each branch in the augmented tree structure takes exactly one character state. Under this process, branches of the same character state share the same OU parameters (sigma). These character states can correspond to the states of an observed discrete character.

This model is equivalent to the MuSSCRat model (May and Moore 2020) if (1) the continuous trait is univariate, and (2) we assume no background rate variation. The major difference between these two models lies in the probability computation. The MuSSCRat model transforms the branch lengths according to the `sigma` parameter and computes the probability using a state-independent algorithm. The phyloBMSD model does not transform the phylogeny and computes the probability using a state-dependent algorithm.

Applications of this model include:

1. State-dependent evolutionary continuous trait evolution conditionally on discrete character history
In this application, the discrete character history (in `simmap` format) is read (using `readCharacterHistory`) and specified in the `characterHistory` argument. The `sigma` arugments represents the state-dependent diffusion parameter that controls continuous trait evolution in different discrete character states. The continuous trait observations are fixed to the tips using `.clamp()`.

2. Joint inference of discrete character history and state-dependent evolutionary continuous trait evolution
In this application, the discrete character evolution is modeled using the distribution `dnPhyloCTMCDASiteIID` (see Landis et al. 2013, May and Moore 2020). The character history is retrieved from the distribution (using `.characterHistories`) and specified in the `characterHistory` argument. Other parameter specifications are same as above.

## authors
Priscilla Lau
## see_also
dnPhyloOUSD
dnPhyloCTMCDASiteIID
readCharacterHistory
## example
    # setup for a two-state phyloOUSD model
    num_states = 2          # 0 and 1 are the only states

    # option 1: conditioning on a fixed character history
    char_hist = readCharacterHistory( simmap_path )[1]
    
    # option 2: joint inference of character history
    Q <- fnJC(num_disc_states)
    X ~ dnPhyloCTMCDASiteIID(tree, Q, branchRates=1, type="Standard")
    
    # set state-dependent OU parameters
    for (i in 1:num_states){
      sigma2[i] ~ dnLognormal(1, 0.587405)
    }
    
    # basic use of the function
    Y ~ dnPhyloOUSD(char_hist, sigma=sigma2^0.5)



## references
- citation: Landis MJ, Matzke NJ, Moore BR, Huelsenbeck JP (2013). Bayesian analysis of biogeography when the number of areas is large. Systematic biology, 62(6):789-804.
  doi: 10.1093/sysbio/syt040
  url: https://academic.oup.com/sysbio/article-abstract/62/6/789/1708738
- citation: May MR, Moore BR (2020). A Bayesian approach for inferring the impact of a discrete character on rates of continuous-character evolution in the presence of background-rate variation. Systematic biology, 69(3):530-544.
  doi: 10.1093/sysbio/syz069
  url: https://academic.oup.com/sysbio/article/69/3/530/5609130
