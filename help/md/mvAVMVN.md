## name
mvAVMVN
## title
## description
The adaptive variance multivariate-normal proposal of Baele et al. 2017, uses MCMC samples to fit covariance matrix to parameters.

After user-defined waiting time, proposes using covariance matrix epsilon * I + (1 - epsilon) * sigmaSquared * empirical_matrix.

Internally transforms variables based on whether variables are (finitely) bounded, strictly positive, or simplexed.

Non-simplex-valued vector random variables are untransformed.

Add random variables to the move directly (e.g. branch_rates[1], not branch_rates).WARNING: Disabling tuning disables both tuning of proposal variance and learning of empirical covariance matrix.

## details
## authors
## see_also
## example
## references
