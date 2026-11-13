# Reader benchmarks

Run the same script sequentially against main and the working branch, using
the same Julia runtime, resolved dependency versions, fixtures and thread count:

```sh
julia --project=/path/to/main benchmark/readers.jl /path/to/shared/test/assets 10
julia --project=. benchmark/readers.jl /path/to/shared/test/assets 10
```

The project directories must contain Project.toml and resolved dependencies.
The script reads all six original fixtures, measuring mesh-only and all-field
reads. It also measures the four original Eclipps surface extraction workloads.
Each case is warmed up before measuring elapsed time, allocated bytes and
allocation counts. Garbage collection runs before each sample, outside the
timed region. Results describe warm filesystem caches, not cold disk throughput
or compilation. Unsupported cases report their errors on stderr and a CSV row.
Cell/node counts help identify changes in the compared workload; successful
benchmarks do not replace correctness tests.

The largest original fixture has 235,008 cells and 296,320 nodes; ten samples
of its full reads take several minutes.

This optional script is outside normal tests and performs no downloads. Large
CGNS example downloads remain future work and would require explicit opt-in.
Record commit IDs and environment alongside any published comparison.
