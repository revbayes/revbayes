## name
mvAdaptiveRJSwitch
## title
Adaptive Reversible-Jump (RJ) move
## description
A move that performs a reversible-jump between a fixed value and a value drawn from a distribution. The standard approach is that the value is drawn from the prior distribution. As this can be inefficient for broad priors, we learn here the proposal distribution. Currently, we only support a normal distribution. Thus, we learn the mean and variance of this normal distribution during the learning phase before actually applying it.
## details
## authors
Sebastian Höhna
## see_also
mvRJSwitch
## example

theta ~ dnReversibleJumpMixture(0.01,
                        dnUnif( 0.0, 0.1 ),
                        0.5)

moves.append( mvAdaptiveRJSwitch(theta, waitBeforeLearning=100,
                                        waitBeforeUsing=1000,
                                        updateEvery=10,
                                        weight=10.0) )
## references
