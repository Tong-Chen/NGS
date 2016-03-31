#!/usr/bin/env python
# Script to select the best alignments in a sorted sam file
# also adds the XS:A attribute if not present returning a SAM file with the 
# necessary 9 cols and the xs attribute (add more?)

import sys, re

def getstrand(bit):
    if bit == '0':
        return 'XS:A:-'
    else:
        return 'XS:A:+'

infile = open(sys.argv[1],'r')

tmp = infile.next().split()
while True:
	try:
		tmp2 = infile.next().split()
	except:
		break
	
	if tmp[0] != tmp2[0]:
		strand = getstrand(tmp[1])
		print '\t'.join(tmp[0:11])+'\t'+strand
		tmp = tmp2
		continue
	else:	
	    
		if 'H' in tmp[5]:
			strand=getstrand(tmp2[1])
			print '\t'.join(tmp2[0:11])+'\t'+strand
			try:				
				tmp = infile.next().split()
			except:
				break

		elif 'H' in tmp2[5]:
			strand=getstrand(tmp[1])
			print '\t'.join(tmp[0:11])+'\t'+strand
			try:
			    tmp = infile.next().split()
			except:
			    break

		elif tmp[4] == '255' and tmp2[4] != '255':
			strand = getstrand(tmp2[1])
			print '\t'.join(tmp2[0:11])+'\t'+strand
			try:				
			    tmp = infile.next().split()
			except:
				break
		
		elif tmp[4] != '255' and tmp2[4] == '255':
			strand = getstrand(tmp[1])
			print '\t'.join(tmp[0:11])+'\t'+strand
			try:
			    tmp = infile.next().split()
			except:
			    break

		elif tmp[4]==tmp2[4]:
			print '\t'.join(tmp[0:11])
			try:
			    tmp = infile.next().split()
			except:
			    break
infile.close()
