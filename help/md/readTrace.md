## name
readTrace
## title
Read an MCMC log file.
## description
Reads parameter values from a log file, usually produced as the output of an MCMC or MCMCMC run.
## details

Read an MCMC log file with field delimited by `separator`.
Then drop the first `burnin` iterations if `burnin` is an integer,
or the fraction `burnin` of iterations if `burnin` if a Real number.
Then we keep every `n`th entry if the `thinning` is `n`.

## authors
## see_also
readAncestralStateTrace
readTreeTrace
readCharacterData
Trace
## example
    # Read in a log file as a vector of traces
    traces = readTrace(filename, burnin = burnin)

## references
