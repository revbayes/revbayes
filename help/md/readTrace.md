## name
readTrace
## title
readTrace
## description
Read an MCMC log file.
## details
Read an MCMC log file with field delimited by `separator`.
Then drop the first `burnin` iterations if `burnin` is an integer,
or the fraction `burnin` of iterations if `burnin` if a Real number.
Then we keep every `n`th entry if the `thinning` is `n`.
## authors
## see_also
## example
trace <- readTrace(filename, burnin=burnin)[1]
## references
