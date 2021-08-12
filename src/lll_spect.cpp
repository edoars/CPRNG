/*  Written in 2019 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* Prints the figures of merit for a multiplier and a power-of-two
   modulus. If MULT is defined, computes figures of merit for an MCG;
   otherwise, for an LCG.

   See also Karl Entacher & Thomas Schell's code associated with the paper

   Karl Entacher, Thomas Schell, and Andreas Uhl. Efficient lattice assessment
   for LCG and GLP parameter searches. Mathematics of Computation,
   71(239):1231–1242, 2002.

   at

   https://web.archive.org/web/20181128022136/http://random.mat.sbg.ac.at/results/karl/spectraltest/
*/

#include <iostream>
#include <limits>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

using namespace NTL;
using namespace std;

#include "common.cpp"

int main(int argc, char *argv[]) {
	if (argc != 5) {
		cerr << "USAGE: " << argv[0] << " LAG MAXDIM MULTIPLIER MODULUS" << endl << endl;
		cerr << "For the given multiplier and power-of-two modulus uses the" << endl;
		cerr << "LLL lattice-reduction algorithm to compute figures of merit" << endl;
		cerr << "up to the specified maximum dimension for the given lag." << endl;
		cerr << "A lag of one gives the standard spectral test. Prints the minimum" << endl;
		cerr << "spectral score, the harmonic spectral score, the multiplier" << endl;
		cerr << "in decimal and hexadecimal, the lag and the figures of merit" << endl;
		cerr << "up to the specified maximum dimension." << endl;
		exit(1);
	}

	const int lag = atoi(argv[1]);

	if (lag < 1) {
		cerr << "The lag must be strictly positive" << endl;
		exit(1);
	}

	const int max_dim = atoi(argv[2]);

	if (max_dim > dim_max) {
		cerr << "Maximum possible number of dimensions: " << dim_max << endl;
		exit(1);
	}

	if (max_dim < 2) {
		cerr << "Minimum possible number of dimensions: 2" << endl;
		exit(1);
	}

	ZZ a = strtoZZ(argv[3]);
	ZZ mod = strtoZZ(argv[4]);

	if (a >= mod) {
		cerr << "The multiplier must be smaller than the modulus" << endl;
		exit(1);
	}

	// Entacher's characterization of lagged lattices (https://dl.acm.org/doi/10.1145/301677.301682)
	ZZ alag = PowerMod(a, lag, mod);
	mod /= GCD(mod, conv<ZZ>(lag));
#ifdef MULT
	mod /= 4;
#endif
	// Compute the normalization factor starting from gamma_t
	for(int d = 2; d <= dim_max; d++)
		norm[d - 2] = to_double(to_RR(1) / (pow(to_RR(norm[d - 2]), to_RR(1./2)) * pow(to_RR(mod), to_RR(1) / to_RR(d))));

	mat_ZZ mat;
	mat.SetDims(max_dim, max_dim);

	double harm_norm = 0, min_fm = numeric_limits<double>::infinity(), harm_score = 0, cur_fm[dim_max];

	for (int d = 2; d <= max_dim; d++) {
		mat.SetDims(d, d);
		// Dual lattice (see Knuth TAoCP Vol. 2, 3.3.4/B*).
		mat[0][0] = mod;
		for (int i = 1; i < d; i++) mat[i][i] = 1;
		for (int i = 1; i < d; i++) mat[i][0] = -power(alag, i);
		ZZ det2;
		// LLL reduction with delta = 0.999999999
		LLL(det2, mat, 999999999, 1000000000);

		double min2 = numeric_limits<double>::infinity();

		for (int i = 0; i < d; i++) {
			//cout << mat[i] << endl;
			min2 = min(min2, to_double(mat[i] * mat[i]));
		}
		cur_fm[d - 2] = norm[d - 2] * sqrt(min2);
		min_fm = min(min_fm, cur_fm[d - 2]);
		harm_score += cur_fm[d - 2] / (d - 1);
		harm_norm += 1. / (d - 1);
	}

	printf("%8.6f\t%8.6f\t", min_fm, harm_score / harm_norm);
	cout << a << "\t" << "0x" << hex(a) << "\t" << lag;
	for (int d = 2; d <= max_dim; d++) printf("\t%8.6f", cur_fm[d - 2]);
	cout << endl;
	return 0;
}
