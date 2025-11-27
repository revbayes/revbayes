## name
readTrace
## title
Read an MCMC log file.
## description
Reads parameter values from a log file, usually produced as the output of an
MCMC or MCMCMC run.
## details
Either a file name or a directory must be provided as argument. If a folder is
provided, all the files that contain trees in that directory are read in
the same object.

When reading individual log files, field is delimited by `separator`. We drop
the first `burnin` iterations if `burnin` is an integer greater than or equal
to 1, or the fraction `burnin` of iterations if `burnin` is a real number
smaller than 1. Then we keep every `n`th entry if the `thinning` is `n`.

`trace = readTrace(..., nruns = 1)` returns a `Trace[]` object.
`trace = readTrace(..., nruns > 1)` returns a `Trace[][]` object, in which
`trace[n]` corresponds to the trace of run `n`.

## authors
## see_also
readAncestralStateTrace
readTreeTrace
Trace
## example
    # read in a single log file as a vector of traces
    traces = readTrace("out.log", burnin=0.5)
    
    # read in two log files called 'out_run_1.log', 'out_run_2.log'
    traces = readTrace("out.log", burnin=0.5, nruns=2)

## references
