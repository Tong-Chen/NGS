#!/usr/bin/env python

#Script to count the number of the four different events characterized by
#astalavista, based on the summary output from summary.py

import sys

infile = open(sys.argv[1],'r').readlines()
exon = 0
acc = 0
don = 0
intron = 0


for line in infile:
	tmp = line.split()
	count = int(tmp[1])
	tmp2 = str(tmp[0]).split(',')
	iso1 = str(tmp2[0])
	iso2 = str(tmp2[1])
	if iso1 == '0' and iso2.endswith('^'):
		exon += count
		#f1 = iso2.replace('^',',')
		#f2 = f1.replace('-',',')
		#length = len(f2.split(','))
		#print length
	elif iso1.endswith('^') and iso2.endswith('^'):
		don += count
	elif iso1.endswith('-') and iso2.endswith('-'):
		acc += count
	elif iso1 == '0' and iso2.endswith('-'):
		intron += count
	else:
		print >>sys.stderr, 'no type'

print 'Exon skipping events: %d' %(exon)
print 'Alternative acceptor sites: %d' %(acc)
print 'Alternative donor sites: %d' %(don)
print 'Intron retention events: %d' %(intron)
print 'All events: %d' %(exon+acc+don+intron)
	
