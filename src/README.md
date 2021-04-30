Code for the Spectral Test
==========================

The code in this directory deals with (approximate) spectral figures of
merit and aggregate scores for LCGs and MCGs with modulus that is a power
of 2. It can be easily adapted to other moduli. It uses
[NTL](https://www.shoup.net/ntl/).

- `lll_search.cpp` code searches for good multipliers for LCGs and MCGs.

- `lll_spect.cpp` prints spectral scores and figures of merit for a given multiplier.

- `lll_printdat.cpp` prints configuration files for
  [LatticeTester](https://github.com/umontreal-simul/latticetester) for a
  given multiplier.

- `benchmark.c` is a simple microbenchmark comparing different multiplier
  sizes (compilation instruction are at the start of the file).

The code for spectral figures of merit is based on the
Lenstra–Lenstra–Lovász lattice basis reduction algorithm, and as such it
computes approximate results. Up to dimension 8, however, the
approximation is usually excellent (see
<https://doi.org/10.1090/S0025-5718-01-01415-6>). If you want to get exact
results for a multiplier, use LatticeTester with the configuration files
provided by `lll_printdat`.

See the `comp.sh` script for the options needed to compile for the 
case of LCGs or MCGs.
