#!/usr/bin/env python

#Script to extract GTF annotation over a specified cuf off value 
#Usage: extract_gtf.py <FPKM-values file> [Cut off] <GTF-file>

import sys, re

infile = open(sys.argv[1],'r')
cutoff = float(sys.argv[2])
gtf = open(sys.argv[3],'r')

isoforms = []

reg = r'transcript_id "(TCONS_[0-9]+)"'

for x in infile:
	try:
		parts = x.split('\t')
		if float(parts[9]) > cutoff or float(parts[13]) > cutoff:
			isoforms.append(parts[0])			
			
	except:
		continue

for y in gtf:
	match = re.search(reg, y)
	trans = match.groups(1)
	if trans[0] in isoforms:
		sys.stdout.write(y)
	else:
		continue


