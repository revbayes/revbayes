### Current changes (development branch)

#### New models/analyses

#### New features

#### Bug fixes


### Version 1.1.1

**Warning**: this version includes changes to the Random Number Generator, meaning analysis output will be different from v1.1.0, even when run with the same seed. 

#### New models/analyses

 * new PoMo rate matrices with selection

#### Bug fixes

 * fixed an error in the likelihood calculation for discrete characters when using ambiguous characters and partitioning standard data by state number
 * fixed a bug causing unbounded likelihoods when using PhyloBrownianProcessREML with unrooted trees
 * fixed a crash when using ordered rate matrices
 * fixed a crash when printing the value of a rate matrix
 * fixed a crash when using nested clade constraints
 * fixed a crash when using DEC models
 * fixed an issue with using simplex indexes in deterministic nodes
 * fixed a bug allowing NodeTimeSlideUniform to move sampled ancestor nodes
