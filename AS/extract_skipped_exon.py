#!/usr/bin/env python

#Script to extract regions containing a skipped exon
#Input file is the raw output from ASTA

import re, sys

data = open(sys.argv[1], 'r')

#structure "0,1-2^"
reg1 = r'structure ("0,1-2\^)"'
#splice_chain ",11691329-11691389^
reg2 = r'splice_chain ",([0-9]+)-([0-9]+)'
#transcript_id "TCONS_00000024
reg3 = r'transcript_id "(TCONS_[0-9]+)'
for x in data:
	if re.search(reg1, x):
		tmp = re.search(reg2, x)
		start = tmp.groups()[0]
		stop = tmp.groups()[1]
		tmp2 = x.split('\t')
		chrom = tmp2[0]
		strand = tmp2[6]
		trans = re.search(reg3,x)
		name = trans.groups()[0]
		
		if strand == '-':
			out = [chrom, stop, start, name,'score',strand,tmp2[3],tmp2[4]]
		else:
			out = [chrom, start, stop, name,'score',strand,tmp2[3],tmp2[4]]
	
		print '\t'.join(out)	
	else:
		continue
		
	
	
