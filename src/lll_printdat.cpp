/*  Written in 2019 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* Writes input files for LatticeTester for a multiplier and a power-of-two
   modulus. If MULT is defined, writes files for an MCG; otherwise, for an LCG. */

#include <iostream>
#include <fstream>
#include <limits>
#include <cstring>
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
		cerr << "USAGE: " << argv[0] << " MAXDIM MULTIPLIER MODULUS BASENAME" << endl << endl;
		cerr << "For the given multiplier and power-of-two modulus prints" << endl;
		cerr << "input files for LatticeTester to compute the spectral" << endl;
		cerr << "figures of merit up to the specified maximum dimension." << endl;
		exit(1);
	}

	const int max_dim = atoi(argv[1]);

	if (max_dim > dim_max) {
		cerr << "Maximum possible number of dimensions: " << dim_max << endl;
		exit(1);
	}

	if (max_dim < 2) {
		cerr << "Minimum possible number of dimensions: 2" << endl;
		exit(1);
	}

	ZZ a = strtoZZ(argv[2]);
	ZZ mod = strtoZZ(argv[3]);
	string basename(argv[4]);

	if (a >= mod) {
		cerr << "The multiplier must be smaller than the modulus" << endl;
		exit(1);
	}

	// Compute the normalization factor starting from gamma_t
	for(int d = 2; d <= dim_max; d++)
		norm[d - 2] = conv<double>(conv<RR>(1) / (pow(conv<RR>(norm[d - 2]), conv<RR>(1./2)) * pow(conv<RR>(mod), conv<RR>(1) / conv<RR>(d))));

	mat_ZZ x;
	x.SetDims(max_dim, max_dim);

#ifdef MULT
	mod /= 4;
#endif

	for (int d = 2; d <= max_dim; d++) {
		ofstream s;
		s.open(basename + "-" + to_string(d) + ".dat");
//		s << "SPECTRAL BESTLAT" << endl << "PREREDDIETER" << endl << d << endl;
		s << "SPECTRAL BESTLAT" << endl << "BKZ ARBITRARY 0.9999999999 100" << endl << d << endl;
		x.SetDims(d, d);
		// Dual lattice (see Knuth TAoCP Vol. 2, 3.3.4/B*).
		x[0][0] = mod;
		for (int i = 1; i < d; i++) x[i][i] = 1;
		for (int i = 1; i < d; i++) x[i][0] = -power(a, i);

		for(int i = 0; i < d; i++) {
			for(int j = 0; j < d; j++) {
				s << x[i][j] << " ";
			}
			s << endl;
		}

		// Max nodes, output
		s << 1000000000 << endl << "TERMINAL" << endl;
		s.close();
	}
}
