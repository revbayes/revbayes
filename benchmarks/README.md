# Benchmarks

## Overview

These benchmarks will depend on whether `useScaling` is set to `true` or `false`

## bench-unrooted — single chain baseline

```
time rb bench-unrooted.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 100000
time rb bench-unrooted.Rev ../tests/data/EF-Tu-12.fasta AA 10000
```

## bench-replicates — Phase 1 OpenMP (parallel replicates)

Measures speedup from parallelising independent replicate analyses.
Run serial baseline (nruns=1) then parallel (nruns=N) and compare wall time.

```
# Serial baseline
time rb bench-replicates.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 10000 1

# 2 / 4 / 8 replicates in parallel (requires -omp true build)
time rb bench-replicates.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 10000 2
time rb bench-replicates.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 10000 4
time rb bench-replicates.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 10000 8
```

Expected: near-linear speedup (0.7–1.0× per added replicate) since runs are independent.

## bench-mc3 — Phase 2 OpenMP (parallel MC³ chains)

Measures speedup from parallelising chains within a single MCMCMC analysis.
Compare nchains=1 vs nchains=N.

```
# Serial baseline (1 chain = plain MCMC)
time rb bench-mc3.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 10000 1

# 2 / 4 / 8 chains in parallel (requires Phase 2 OpenMP + -omp true build)
time rb bench-mc3.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 10000 2
time rb bench-mc3.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 10000 4
time rb bench-mc3.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 10000 8
```

Expected: near-linear speedup once Phase 2 (Mcmcmc.cpp chain loop) is parallelised.
