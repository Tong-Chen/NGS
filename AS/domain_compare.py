#!/usr/bin/env python

#Script to compare the functional domains of known and novel isoforms
#Input; 1- File with novel-known pairs, 2- domains novel, 3- domains known transcripts

import sys

pair = open(sys.argv[1],'r').readlines()
novel = open(sys.argv[2],'r').readlines()
known = open(sys.argv[3],'r').readlines()

def complists(novellist, knownlist):
	retained = []
	gain_off = []
	loss_off = []

	for i in novellist:
		if i in knownlist:
			retained.append(i)
		elif i not in knownlist:
			gain_off.append(i)
	for z in knownlist:
		if z not in novellist:
			loss_off.append(z)
		else:
			continue
	return retained, gain_off, loss_off
lost_all = 0
lost_some = 0
retained_all = 0
retained_some = 0
gained = 0
other = 0
none = 0
	
for line in pair:
	ndom = []
	kdom = []
	k, n = line.split()
	for x in novel:
		if n in x:
			tmp = x.split('\t')
			ndom.append(tmp[6])
		else:
			continue
		
	for y in known:
		if k in y:
			tmp2 = y.split('\t')
			kdom.append(tmp2[7])

	out1, out2, out3 = complists(ndom, kdom)

	if len(kdom) == 0 and len(ndom) == 0:
		none += 1
	elif len(kdom) > 0 and len(ndom) == 0:
		lost_all += 1
	elif len(ndom) > 0 and len(out3) > 0:
		lost_some += 1
	elif len(out1) > 0 and len(out3) == 0:
		retained_all += 1
	elif len(out2) > 0:
		gained += 1
		print n,'\t',k,'---',','.join(out1),'---',','.join(out2),'---',','.join(out3)
	else:
		other += 1

#print 'none: %d, lost all: %d, lost_some: %d, retained_all: %d, gained: %d, other: %d' %(none, lost_all, lost_some, retained_all, gained, other)

	#print n,'\t',k,'---',','.join(out1),'---',','.join(out2),'---',','.join(out3)
	
