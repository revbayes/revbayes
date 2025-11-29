## name
dnPhyloOrnsteinUhlenbeckStateDependent (dnPhyloOUSD)
## title
Phylogenetic state-dependent Ornstein-Uhlenbeck process
## description
Univariate Ornstein-Uhlenbeck process over an augmented phylogeny (i.e., discrete character history)
## details
The phylogenetic state-dependent Ornstein-Uhlenbeck process is an extension toto Hansen's (1997) Ornstein-Uhlenbeck process with state-dependent optima.
The probability is computed using a pruing algorithm, which is derived from FitzJohn's (2012) pruning algorithm for a phylogenetic state-independent Ornstein-Uhlenbeck process.
Specifically, `dnPhyloOUSD` uses an augmented tree structure, i.e., the character history with different character states (also referred to as "regimes"; Hansen 1997).
Shifts in the character state is discrete and does not occur concurrently with speciation events, and can be represented by nodes of degree 2.
As such, each branch in the augmented tree structure takes exactly one character state.
Under this process, branches of the same character state share the same OU parameters (alpha, theta, sigma).
These character states can correspond to the states of an observed discrete character.

Applications of this model include:

1. State-dependent evolutionary continuous trait evolution conditionally on discrete character history

In this application, the discrete character history (in `simmap` format) is read (using `readCharacterHistory`) and specified in the `characterHistory` argument.
The `alpha`, `theta`, and `sigma` arugments represent the state-dependent OU parameters that controls continuous trait evolution in different discrete character states.
Furthermore, the assumption of the continuous trait value at the root is specified using the `rootTreatment` argument.
The continuous trait observations are fixed to the tips using `.clamp()`.

2. Joint inference of discrete character history and state-dependent evolutionary continuous trait evolution
In this application, the discrete character evolution is modeled using the distribution `dnPhyloCTMCDASiteIID` (see Landis et al. 2013, May and Moore 2020).
The character history is retrieved from the distribution (using `.characterHistories`) and specified in the `characterHistory` argument.
The `alpha`, `theta`, and `sigma` arugments represent the state-dependent OU parameters that controls continuous trait evolution in different discrete character states.
Furthermore, the assumption of the continuous trait value at the root is specified using the `rootTreatment` argument.
The continuous trait observations are fixed to the tips using `.clamp()`.

## authors
Priscilla Lau
## see_also
dnPhyloBMSD
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
      theta[i] ~ dnUniform(-10, 10)
      alpha[i] ~ dnLognormal(ln(2), 0.587405)
      sigma2[i] ~ dnLognormal(1, 0.587405)
    }
    
    # basic use of the function (assuming the continuous trait value at the root is the same at the state-specific optimum)
    Y ~ dnPhyloOUSD(char_hist, theta=theta, rootTreatment="optimum", alpha=alpha, sigma=sigma2^0.5)



## references
- citation: Hansen TF (1997). Stabilizing selection and the comparative analysis of adaptation. Evolution, 51(5):1341-1351.
  doi: 10.1111/j.1558-5646.1997.tb01457.x
  url: https://academic.oup.com/evolut/article/51/5/1341/6757302
- citation: FitzJohn RG (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution, 3(6):1084-1092.
  doi: 10.1111/j.2041-210X.2012.00234.x
  url: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2012.00234.x
- citation: Landis MJ, Matzke NJ, Moore BR, Huelsenbeck JP (2013). Bayesian analysis of biogeography when the number of areas is large. Systematic biology, 62(6):789-804.
  doi: 10.1093/sysbio/syt040
  url: https://academic.oup.com/sysbio/article-abstract/62/6/789/1708738
- citation: May MR, Moore BR (2020). A Bayesian approach for inferring the impact of a discrete character on rates of continuous-character evolution in the presence of background-rate variation. Systematic biology, 69(3):530-544.
  doi: 10.1093/sysbio/syz069
  url: https://academic.oup.com/sysbio/article/69/3/530/5609130
