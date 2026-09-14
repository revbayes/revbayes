# Benchmark results

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

##### AVX  (-march=native)  (1.32x slower)
2.090
2.093
 2.093
2.102
2.111
