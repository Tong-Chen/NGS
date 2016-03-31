#!/usr/bin/env python
# Arguments: 1: name of exon list (in bed format), 2: number of exons to include, 3: Bedfile, 4: Number of permutations

import sys, subprocess

usage = 'Usage: command exon_list(bedfile) number_of_exons_per_permutation(int)	file_of_peaks(bedfile) number_of_permutations'

#Location of bedtools
loc = '../../programs/BEDTools-Version-2.15.0/bin/bedtools'

try:
	name = sys.argv[1]
	number = str(sys.argv[2])
	bedfile = sys.argv[3]
	perm = sys.argv[4]
	count = 0
except:
	print usage
	sys.exit(1)


for n in range(int(perm)):
	tmp = subprocess.Popen(['shuf', '-n',number, name], stdout = subprocess.PIPE)

	u = subprocess.Popen([loc, 'intersect', '-a', 'stdin', '-b', bedfile], stdin = tmp.stdout, stdout = subprocess.PIPE)

	t = subprocess.Popen(['wc','-l'],stdin = u.stdout, stdout = subprocess.PIPE)

	count += int(t.communicate()[0])



out = float(count)/float(perm)
print out

