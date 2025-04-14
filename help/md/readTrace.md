## name
readTrace
## title
Read a trace file, usually produced as the output of an MCMC or MCMCMC run.
## description
Reads parameter values from a file, perhaps containing the results of an MCMC run, and saves them in one object.
## details
Given a delimiter-separated table of values, `readTrace` stores each column as an entry of a `Trace[]` vector, such that `Trace[1]` contains every `thinning`-th value from the first column of `file`, after excluding the first `burnin` values.

Pertinent methods of a `Trace` object include:

- `trace.getBurnin()`: Return the number of samples that were discarded as burnin, after thinning.
- `trace.getValues()`: Return a vector (e.g. `Real[]`) containing the values of the sampled parameter,
  after excluding burnin samples.
- `trace.setBurnin()`: Intended to allow modification of the burnin size. Currently has no effect.
- `trace.size()`, `trace.getNumberSamples()`: Report the number of values stored in `trace`,
  including (`post = FALSE`) or excluding (`post = TRUE`) burnin samples.
- `trace.summarize()`: Display summary statistics of trace.

## authors
## see_also
readTreeTrace
readCharacterData
## example
## references
