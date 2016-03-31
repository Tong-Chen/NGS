#!/usr/bin/env python

#Script to get all TSS's in from the GFF file (GTF2GFF using GFFread) from cufflinks.
#Outputs unique TSS

import sys

import sys,re

data = open(sys.argv[1],'r')
tss = {}
reg = r'geneID=(XLOC_[0-9]+)'



for x in data:
        tmp = x.split()
        if 'transcript' in x:
                if tmp[2] == 'transcript':
                        name = re.search(reg, x)
                        name2=name.groups()[0]
                        if tmp[6]=='+':
                                name3 = name2+'*'+tmp[3]
				start = tmp[3]
				stop = int(tmp[3])+1
                        else:
                                name3 = name2+'*'+tmp[4]
				start = tmp[4]
                                stop = int(tmp[4])+1

                        if name3 in tss:
                                continue
                        else:
				#print BED-info, name == xloc_123*start
				out = ['chr'+str(tmp[0]),str(start),str(stop),name3,'score',tmp[6]]
				print '\t'.join(out)

#				print 'chr'+tmp[0],'\t',start,'\t',stop,'\t',name3,'\t','score','\t',tmp[6]
			        tss[name3]=1
				

        else:
                continue




