Spectral test search
====================

The code in this directory computes approximate spectral scores for LCGs
and MLCGs with modulus that is a power of 2. It can be easily adapted to
other moduli.

- `lll_search.cpp` code searches for good multipliers for LCGs and MLCGs.
- `lll_spect.cpp` prints spectral scores for a given multiplier.
- `lll_printdat.cpp` prints configuration files for (LatticeTester)[https://github.com/umontreal- simul/latticetester] for a given multiplier.

The code is based on the Lenstra–Lenstra–Lovász lattice basis reduction
algorithm, and as such it computes approximate results. Up to dimension 8,
however, the approximation is be excellent (see
<https://doi.org/10.1090/S0025-5718-01-01415-6>). 
It can be easily adapted to other moduli.

See the `comp.sh` script for the options needed to compile for the 
case of LCGs or MCGs.
