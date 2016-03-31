#!/usr/bin/env python

#Script to retrieve the last exon of a locus

import sys,re

data = open(sys.argv[1],'r')
reg = r'(XLOC_[0-9]+)'

for n in data:
	if 'transcript' in n or 'exon' in n:
                parts = n.split()
                if parts[2]== 'transcript':
			n = re.search(reg, n)
                        name = n.groups()
			if parts[6] == '+':
	                        end = parts[4]
			else:
				end = parts[3]
                        lock = 1

		elif parts[6] == '+' and parts[2] == 'exon' and parts[4] == end and lock == 1:
        	                start = parts[3]
                	        lock = 0
                        	print 'chr'+parts[0],'\t',start,'\t', end,'\t', name[0],'\t', 'score', '\t',parts[6]

		elif parts[6] == '-' and parts[2] == 'exon' and parts[3] == end and lock == 1:
				start = parts[4]
				lock = 0
				print 'chr'+parts[0],'\t', end,'\t', start,'\t', name[0],'\t','score','\t', parts[6]
			
		else:
			continue
		
	else:
		continue
