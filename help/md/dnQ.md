## name
dnQ

## title
Distribution on reversible and non-reversible nucleotide rate matrices

## description
Defines a nucleotide substitution rate matrix that can move between
time-reversible and non-reversible models during reversible-jump MCMC.

## details
Use this distribution with `mvMPQRateMatrix` to sample the rate matrix
and propose changes between the two models.

`alpha` is a vector of four positive Dirichlet concentration parameters
for the stationary nucleotide frequencies. Setting all four entries to
one gives a uniform prior on these frequencies.

`fixedModel` defaults to `"free"`, allowing transitions between models.
Set it to `"reversible"` or `"nonreversible"` to restrict the analysis
to that model; model-switching proposals are then disabled.

`lnPriorOdds` specifies log(rho_R/rho_N), where rho_R and rho_N are the
prior probabilities of the reversible and non-reversible models.
Its default, zero, gives equal prior probabilities. Positive values
favor reversibility; negative values favor non-reversibility.

For a free-model analysis with fixed prior odds, the Bayes factor
favoring the non-reversible model is

    BF_NR = (p_N / p_R) * exp(lnPriorOdds)

where p_N and p_R are the posterior model probabilities estimated from
the sampled model indicators. Use production samples after burn-in,
with adequate mixing between models.

If `mvMPQRateMatrix` adapts the model prior during burn-in, freeze that
adaptation before collecting production samples and use the final
value reported by `Q.lnPriorOdds()` in the calculation.

## authors
John P. Huelsenbeck
Sebastian Höhna

## see_also
mvMPQRateMatrix
mvMPQTree
mvReversibilityAwareBranchLengthScale
dnPhyloCTMC

## example
    # Create a rate matrix with equal prior probabilities for the models.
    Q ~ dnQ(alpha=v(1,1,1,1), lnPriorOdds=0.0)

    moves = VectorMoves()
    moves.append(mvMPQRateMatrix(Q, tune=FALSE, tuneModelPrior=FALSE))

    # Define quantities that can be monitored during MCMC.
    reversible := ifelse(Q.isReversible(), 1.0, 0.0)
    frequencies := Q.getPi()
    log_prior_odds := Q.lnPriorOdds()

    # Alternatively, restrict a separate matrix to one model.
    Q_NR ~ dnQ(alpha=v(1,1,1,1), fixedModel="nonreversible")
