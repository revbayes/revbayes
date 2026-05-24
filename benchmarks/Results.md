# Benchmark results

## Phase 3 — site-pattern likelihood parallelisation | sriramv | 2026-05-23

### bench-site-likelihoods.Rev — primates_cytb (23 taxa, 1141 sites, gamma+4-rates, 20k gens)

```
OMP_NUM_THREADS=1   4.87 s   1.00x
OMP_NUM_THREADS=4   9.29 s   0.52x  ← overhead dominates
OMP_NUM_THREADS=8  17.18 s   0.28x  ← overhead dominates
```

**Initial threshold N>64 — regression observed:**
```
OMP_NUM_THREADS=1   4.87 s   1.00x
OMP_NUM_THREADS=4   9.29 s   0.52x  ← fork/join overhead dominates
OMP_NUM_THREADS=8  17.18 s   0.28x  ← fork/join overhead dominates
```

**After raising threshold to N>2048 — regression eliminated:**
```
OMP_NUM_THREADS=1   4.49 s   1.00x
OMP_NUM_THREADS=4   4.24 s   1.06x  ✅ no regression
OMP_NUM_THREADS=8   4.62 s   0.97x  ✅ no regression
```

**Analysis:** Dataset has ~500 compressed unique patterns < 2048 threshold, so the guard
correctly falls through to serial. Fork/join cost (~20 µs per spawn) exceeds per-call
computation at this size. Parallel path activates on large phylogenomic datasets
(≥2048 unique site patterns, many taxa) where work dominates synchronisation overhead.



This is an unofficial place to collect results.
The results are not official.

## Machines

### Ben | centromere | 4/17/2026

#### time rb -o useScaling=false bench-unrooted.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 100000

##### Plain C
1.942
1.956
 1.956
1.961
1.970

##### Plain C / -march=native
1.641
1.641
 1.644
1.645
1.654

##### SSE
2.384
2.382
 2.390
2.394
2.438

##### AVX / -march=native
1.574
1.575
 1.575
1.594
1.596


#### time rb -o useScaling=true bench-unrooted.Rev ../tests/data/primates_and_galeopterus_cytb.nex DNA 100000

##### Plain C (1.52x slower)
2.951
2.958
 2.982
2.983
2.991

##### Plain C / -march=native (1.79x slower)
2.927
2.937
 2.941
2.949
2.960

##### SSE (1.29x slower)
3.066
3.078
 3.081
3.082
3.089

##### AVX  (-march=native)  (1.32x slower)
2.090
2.093
 2.093
2.102
2.111
