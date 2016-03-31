#!/usr/bin/env python

#Script to get the TES that is furthest away (from the TSS) from a locus, from a GFF file (input file == all.GFF)

import sys,re

data = open(sys.argv[1],'r')
tes = {}
reg = r'geneID=(XLOC_[0-9]+)'


#Get all individual TSS's, the dict shows number of transcripts per tes 
for x in data:
	tmp = x.split('\t')
	if 'transcript' in x:
		if tmp[2] == 'transcript':
			name = re.search(reg, x)
			name2=name.groups()[0]
			if tmp[6]=='+':
				name3 = name2+'*'+tmp[4]
			else:
				name3 = name2+'*'+tmp[3]
		
			if name3 in tes:
				tes[name3]+=1
			else:
				tes[name3]=1


	else:
		continue	


#This prints number of transcripts per TES
#for h in tes:
#	print h ,'\t',tes[h]
print 'Number of TES: %d' %(len(tes))

loci = {}

for y in tes:
	name, start = y.split('*')
	if name in loci:
		loci[name]+=1
	else:
		loci[name]=1

#This prints the number of TES per locus
#for i in loci:
#	print i,'\t', loci[i]
print 'Number of loci: %d' %(len(loci))

	
one = 0
two = 0
three = 0
more = 0

for i in loci:
        if loci[i] == 1:
                one += 1
        elif loci[i] == 2:
                two += 1
        elif loci[i] == 3:
                three += 1
        elif loci[i] > 3:
                more += 1

print 'one: %d' %(one)
print 'Two: %d' %(two)
print 'Three: %d' %(three)
print 'More than three: %d' %(more)

