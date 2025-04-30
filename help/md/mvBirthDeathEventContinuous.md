## name
mvBirthDeathEventContinuous
## title
Birth-death proposal for episodic shift in diversification rate
## description
The move mvBirthDeathEventContinuous is designed to propose changes to the number and placement of diversification rate shift
## details
This proposal adds or removes a single continuous event (e.g., diversification rate shift) in a stochastic character map along a phylogenetic tree. At each MCMC iteration, the move randomly choose to either insert a new event (birth) at a random time on a branch or delete an existing event (death). The probability of choosing a birth or death on the current number of events. Each event typically represents a cahnege in continuous parameters, and this move enables reversible-jump MCMC to comare models with different sifts.
## authors
## see_also
## example
## references
