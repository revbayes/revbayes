## name
ConditionalPosteriorOrdinate
## title
Conditional posterior ordinate (a.k.a. cross-validation)
## description
Model selection via leave one out cross-validation. Cross-­ validation assesses the fit of a model by using a subset of the data to estimate parameters (i.e., train the model) and the remaining observations to evaluate the predictive fitness. In the most extreme case of leave-­one-­out cross-­validation, we use all but one observation to estimate the parameters, P(𝜃| X(−i)), then compute the probability of observing the removed observed, P(Xi| 𝜃) , integrate over all parameters, P(Xi| X(−i)) = ∫ p(Xi| 𝜃)p(𝜃| X(−i))d𝜃, and repeat this process for all data points.
## details
A cross-validation analysis assumes that one has a trace of samples of the probabilities for each data point (e.g., site or column) from your data stored in a file. Each data point probability is important for the computation of the leave-one-out cross-validation probability. We read in this trace using the function `ConditionalPosteriorOrdinate`. This "constructor" function requires the `filename` as an argument.

Then, you can calculate the leave-one-out cross-validation probability using the member method `.predictiveProbability()`. In the current implementation, the member method `.predictiveProbability()` requires two argument, the `counts` which are a vector of observations as real numbers, `log` which tells if the probabilities in the trace are log-transformed . The result of this function is the leave-one-out cross-validation probability.
## authors
Sebastian Höhna
## see_also
pathSampler
steppingStoneSampler
## example
obs_sfs = [ 305082, 44248, 32223, 28733, 28220, 26205, 27477, 26618, 27533, 26945, 28736, 28671, 31277, 31250, 34352, 34859, 38331 ]

cpo = ConditionalPosteriorOrdinate( filename="output/StairwayPlot_esfs.log" )
cpo.predictiveProbability( obs_sfs, log=FALSE )
## references
- citation: Lewis PO, Xie W, Chen M-H, Fan Y, Kuo L (2014). Posterior predictive Bayesian phylogenetic model selection. Systematic Biology, 63(3):309-321.
  doi: 10.1093/sysbio/syt068
  url: https://pmc.ncbi.nlm.nih.gov/articles/PMC3985471/pdf/syt068.pdf
- citation: Lartillot N (2023). Identifying the best approximating model in Bayesian phylogenetics: Bayes factors, cross-validation or wAIC? Systematic Biology, 72(3):616-638.
  doi: 10.1093/sysbio/syad004
  url: https://pmc.ncbi.nlm.nih.gov/articles/PMC10276628/pdf/syad004.pdf
- citation: Höhna S, Catalán A (2025). Bayesian StairwayPlot for inferring single population demographic histories from site frequency spectra. Molecular Ecology Resources, 25(6):e14087.
  doi: 10.1111/1755-0998.14087
  url: https://onlinelibrary.wiley.com/doi/pdf/10.1111/1755-0998.14087
