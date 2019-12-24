/*  Written in 2019 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* Prints the figures of merit for a multiplier and a power-of-two
   modulus. If MULT is defined, computes figures of merit for an MCG;
   otherwise, for an LCG. */

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
	if (argc != 4) {
		cerr << "USAGE: " << argv[0] << " MAXDIM MULTIPLIER MODULUS" << endl << endl;
		cerr << "For the given multiplier and power-of-two modulus uses the" << endl;
		cerr << "LLL lattice-reduction algorithm to compute figures of merit" << endl;
		cerr << "up to the specified maximum dimension. Prints the minimum" << endl;
		cerr << "spectral score, the harmonic spectral score, the multiplier" << endl;
		cerr << "in decimal and hexadecimal, the figures of merit up to the" << endl;
		cerr << "specified maximum dimension and lambda (TAB-separated)." << endl;
		exit(1);
	}

	const int max_dim = atoi(argv[1]);

	ZZ a = strtoZZ(argv[2]);
	ZZ mod = strtoZZ(argv[3]);

	if (a >= mod) {
		cerr << "The multiplier must be smaller than the modulus" << endl;
		exit(1);
	}

	if (max_dim > dim_max) {
		cerr << "Maximum possible number of dimensions: " << dim_max << endl;
		exit(1);
	}

	if (max_dim < 2) {
		cerr << "Minimum possible number of dimensions: 2" << endl;
		exit(1);
	}

#ifdef MULT
	mod /= 4;
#endif
	RR sqrt_mod = sqrt(to_RR(mod));

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
		for (int i = 1; i < d; i++) mat[i][0] = -power(a, i);
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

	const double lambda = to_double(to_RR(a) / sqrt_mod);
	printf("%8.6f\t%8.6f\t", min_fm, harm_score / harm_norm);
	cout << a << "\t" << "0x" << hex(a);
	for (int d = 2; d <= max_dim; d++) printf("\t%8.6f", cur_fm[d - 2]);
	printf(lambda < 10000 ? "\t%11.6f\n" : "\t%g\n", lambda); 
	return 0;
}
