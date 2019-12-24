#!/bin/bash

g++ -std=c++17 -O3 -march=native lll_search.cpp -o lll_search -lntl
g++ -std=c++17 -O3 -march=native lll_search.cpp -DMULT -o lll_msearch -lntl
g++ -std=c++17 -O3 -march=native lll_spect.cpp -o lll_spect -lntl
g++ -std=c++17 -O3 -march=native lll_spect.cpp -DMULT -o lll_mspect -lntl
g++ -std=c++17 -O3 -march=native lll_printdat.cpp -o lll_printdat -lntl
g++ -std=c++17 -O3 -march=native lll_printdat.cpp -DMULT -o lll_mprintdat -lntl
