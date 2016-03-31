#!/usr/bin/env python
#Script to count the number of occurences of AS events in a ASTA output file. The output from the script contains the code (see Foissac and Sammeth, #2007) and the count for this type of event 

import sys, re

file = open(sys.argv[1],'r')
reg = r'structure "(\S+)"'

summary = {}


for x in file:
	match = re.search(reg, x)
	out = match.groups(1)
	if out[0] in summary:
		summary[out[0]]+=1
	else:
		summary[out[0]] = 1

for x in summary:
	print x, summary[x]
