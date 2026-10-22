## name
mvMPQRateMatrix

## title
Rate-matrix updates and reversible-jump proposals for dnQ

## description
Updates a nucleotide rate matrix drawn from `dnQ`, including proposals
between the time-reversible and non-reversible models.

## details
`Q` must be a stochastic rate matrix drawn from `dnQ`. When that
distribution uses `fixedModel="free"`, this move can change the model
indicator. When a model is fixed, model-switching proposals are disabled.

`polyhedronAlpha` controls the concentration of the non-reversible
jump proposal. Its default value, one, gives a uniform draw over the
polyhedron. Larger values concentrate the proposal near its center;
values below one favor points nearer its boundary. This parameter
changes the proposal distribution, not the target distribution.

`polyhedronCenter` optionally specifies three coordinates for the
proposal center. An empty vector, the default, uses the reversible
point (0.5, 0.5, 0.5). A center can be chosen using the posterior mean
of `Q.getU()` from a preliminary non-reversible analysis. If a proposed
polyhedron cannot accommodate the requested center, the center is
moved toward the reversible point for that proposal.

`tune` defaults to TRUE and enables tuning when requested by the MCMC.
`tunePolyhedronAlpha`, which defaults to FALSE, enables adaptation of
the proposal concentration during burn-in.

`tuneModelPrior` defaults to FALSE. If enabled together with `tune`,
it adapts the model prior during burn-in to encourage balanced model
occupancy. This changes the target distribution. Finish all adaptation
before collecting production samples, and monitor `Q.lnPriorOdds()`
to retain the final prior odds needed to calculate the Bayes factor.

`weight` controls how frequently this move is selected relative to
other moves.

## authors
John P. Huelsenbeck
Sebastian Höhna

## see_also
dnQ
mvMPQTree
mvReversibilityAwareBranchLengthScale

## example
    Q ~ dnQ(alpha=v(1,1,1,1), lnPriorOdds=0.0)
    moves = VectorMoves()

    # Default proposal with adaptation disabled.
    moves.append(mvMPQRateMatrix(Q, weight=10,
                                tune=FALSE, tuneModelPrior=FALSE))

    # Alternative proposal with a specified center and concentration.
    # In an analysis, choose the proposal settings appropriate to the data.
    shifted_move = mvMPQRateMatrix(Q, weight=10,
                                  tune=FALSE, tuneModelPrior=FALSE,
                                  polyhedronCenter=v(0.56,0.54,0.34),
                                  polyhedronAlpha=10.0)
