#!/usr/bin/python3

#  Written in 2019 by Sebastiano Vigna
#
# To the extent possible under law, the author has dedicated all copyright
# and related and neighboring rights to this software to the public domain
# worldwide. This software is distributed without any warranty.
#
# See <http://creativecommons.org/publicdomain/zero/1.0/>.

import math
import sys
import pandas as pd 

files = [
	["LCG-32-16", 3],
	["LCG-32-17", 3],
	["LCG-32-18", 3],
	["LCG-32-19", 3],
	["LCG-32-24", 3],
	["LCG-32-32", 0],
	["LCG-64-32", 3],
	["LCG-64-33", 3],
	["LCG-64-34", 3],
	["LCG-64-35", 3],
	["LCG-64-48", 0],
	["LCG-64-64", 0],
	["LCG-128-64", 3],
	["LCG-128-65", 3],
	["LCG-128-66", 3],
	["LCG-128-67", 3],
	["LCG-128-68", 3],
	["LCG-128-69", 3],
	["LCG-128-70", 3],
	["LCG-128-71", 3],
	["LCG-128-72", 3],
	["LCG-128-80", 3],
	["LCG-128-96", 0],
	["LCG-128-128", 0],
	["MCG-32-15", 3],
	["MCG-32-16", 3],
	["MCG-32-17", 3],
	["MCG-32-18", 3],
	["MCG-32-19", 3],
	["MCG-32-24", 3],
	["MCG-32-32", 0],
	["MCG-64-31", 3],
	["MCG-64-32", 3],
	["MCG-64-33", 3],
	["MCG-64-34", 3],
	["MCG-64-35", 3],
	["MCG-64-48", 0],
	["MCG-64-64", 0],
	["MCG-128-63", 3],
	["MCG-128-64", 3],
	["MCG-128-65", 3],
	["MCG-128-66", 3],
	["MCG-128-67", 3],
	["MCG-128-68", 3],
	["MCG-128-69", 3],
	["MCG-128-70", 3],
	["MCG-128-71", 3],
	["MCG-128-72", 3],
	["MCG-128-80", 3],
	["MCG-128-96", 0],
	["MCG-128-128", 0],
]

for p in files:

    file = p[0]
    hi_mask = p[1]
    data = pd.read_csv(file, sep='\t', names = ['M₈', 'H₈', 'd', 'h', 'f₂', 'f₃', 'f₄', 'f₅', 'f₆', 'f₇', 'f₈']) 

    if hi_mask != 0:
        mask = hi_mask << (math.floor(math.log(int(data['d'][0]), 2)) - math.floor(math.log(hi_mask, 2)))
        # Manually select lines satisfying the mask (too large integers for pandas)
        loc = []
        for i in range(len(data)):
            if (int(data['d'][i]) & mask) == mask:
                loc.append(i)
        data = data.loc[loc]

    lowerm = data['M₈'].quantile(.999)
    lowerh = data['H₈'].quantile(.999)
    lowerm2 = data['M₈'].quantile(.9999)
    lowerh2 = data['H₈'].quantile(.9999)
    lowerm3 = data['M₈'].quantile(.99999)
    lowerh3 = data['H₈'].quantile(.99999)
    lowerm4 = data['M₈'].quantile(.999999)
    lowerh4 = data['H₈'].quantile(.999999)

    besth = data.sort_values('H₈', ascending=False).iloc[0]['h']
    sys.stdout.write("%s\t%s\n" % (file, besth))
    row = data[data['M₈'] >= lowerm].sort_values('H₈', ascending=False).iloc[0]
    besth2 = row['h']
    sys.stdout.write("%s\t" % file)
    s = ""
    if besth == besth2:
        row = data[data['M₈'] >= lowerm2].sort_values('H₈', ascending=False).iloc[0]
        besth2 = row['h']
        s = "\t${}^*$"
        if besth == besth2:
            row = data[data['M₈'] >= lowerm3].sort_values('H₈', ascending=False).iloc[0]
            besth2 = row['h']
            s = "\t${}^{**}$"
            if besth == besth2:
                row = data[data['M₈'] >= lowerm4].sort_values('H₈', ascending=False).iloc[0]
                besth2 = row['h']
                s = "\t${}^{***}$"

    sys.stdout.write("%s%s\n\n\n" % (besth2, s))

    if len(data[(data['M₈'] >= row['M₈']) & (data['H₈'] >= row['H₈']) & ((data['M₈'] > row['M₈']) | (data['H₈'] > row['H₈']))].index) != 0:
        print("Non-Pareto optimal harmonic multiplier ", besth2)

    bestm = data.sort_values('M₈', ascending=False).iloc[0]['h']
    sys.stderr.write("%s\t%s\n" % (file, bestm))
    row = data[data['H₈'] >= lowerh].sort_values('M₈', ascending=False).iloc[0]
    bestm2 = row['h']
    sys.stderr.write("%s\t" % file)
    s = ""
    if bestm == bestm2:
        row = data[data['H₈'] >= lowerh2].sort_values('M₈', ascending=False).iloc[0]
        bestm2 = row['h']
        s = "\t${}^*$"
        if bestm == bestm2:
            row = data[data['H₈'] >= lowerh3].sort_values('M₈', ascending=False).iloc[0]
            bestm2 = row['h']
            s = "\t${}^{**}$"
            if bestm == bestm2:
                row = data[data['H₈'] >= lowerh4].sort_values('M₈', ascending=False).iloc[0]
                bestm2 = row['h']
                s = "\t${}^{***}$"
    sys.stderr.write("%s%s\n\n\n" % (bestm2, s))

    if len(data[(data['M₈'] >= row['M₈']) & (data['H₈'] >= row['H₈']) & ((data['M₈'] > row['M₈']) | (data['H₈'] > row['H₈']))].index) != 0:
        print("Non-Pareto optimal minimal multiplier ", bestm2)

