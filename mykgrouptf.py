#!/usr/bin/python

import sys
import nwalign as nw

if len(sys.argv) != 2:
    print >>sys.stderr, 'wrong, using %s filename' % sys.argv[0]
    sys.exit(0)

filein=open(sys.argv[1], 'r');
fileout=open(sys.argv[1]+".grouped", 'w');

threshold=0.7
i=0
rep=[]
locus_num=[]
locus=[]
rep_count=[]
rep_pos=[]

curloc='';
for line in filein.readlines():
    if line[0] == '>':
        i=i+1
        curloc=line[1:-1]
    else:
        rep.append(line[:line.find(':')])
        singreps=line.split("#")
        rep_count.append(len(singreps)-1)
        poslist=[]
        for singrep in singreps:
            if singrep != "\n":
                repinfo=singrep.split(":")
                poslist.append(int(repinfo[1]))
        rep_pos.append(poslist)                
        locus_num.append(i)
        locus.append(curloc)
        
fullscore=0
array=()
score=0

for k in range(len(locus)):
    print "Processing repeat", k, "of", locus[k]
    fullscore=nw.score_alignment(rep[k], rep[k], gap_open=-5,\
        gap_extend=-2, matrix='/home/CT/server/pybin/BLOSUM62')
    print >> fileout, ">", k, locus[k], rep[k], rep_count[k], rep_pos[k]
    for j in range(len(rep)):
        array=nw.global_align(rep[k], rep[j], gap_open=-5,\
            gap_extend=-2, matrix='/home/CT/server/pybin/BLOSUM62')
        score=nw.score_alignment(array[0], array[1], gap_open=-5,\
            gap_extend=-2, matrix='/home/CT/server/pybin/BLOSUM62')
        if score>0 and score/float(fullscore)>=threshold and j!=k:
            print >> fileout, j, locus[j], rep[j], rep_count[j], rep_pos[j]

filein.close()
fileout.close()
