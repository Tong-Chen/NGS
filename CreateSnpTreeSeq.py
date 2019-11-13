#!/usr/bin/env python
#coding=utf-8
#这个程序写的很乱，这个程序就是把很多vcf求并集，把snp连成序列，去掉indel，去掉逗号的多突变
#掉逗号的多突变
import re,os,sys,glob
from Bio import SeqIO

if len(sys.argv) < 3:
    print >>sys.stderr, sys.argv[0]+" vcf_list_file vcf_suffix" 

readme=sys.argv[1]
Suffix=sys.argv[2]

#FastaFile=sys.argv[2]
d={}
for i in open(readme):
    i=i.rstrip()
    f=i+Suffix
    i=os.path.basename(i)
    d[i]={}
    #print ">"+".".join(f.split(".")[:-1])
    for x in open(f):
        x=x.rstrip()
        m=re.search("^#",x)
        if m:
            continue
        l=x.split("\t")
        if "INDEL" in l[7]:
            continue
        if len(l[3])!=1 or len(l[4])!=1 or l[3]=="." or l[4]=="." or l[3]=="N" or l[4]=="N":
            continue
        if not "," in l[4]:
            try:
                if not int(l[1]) in d[i][l[0]]:
                    d[i][l[0]][int(l[1])]=(l[3],l[4])
                else:
                    pass
            except:
                d[i][l[0]]={}
                d[i][l[0]][int(l[1])]=(l[3],l[4])
#scaffold={}
position={}
for f in d:
    for s in d[f]:
        #try:
        #    scaffold[s]+=d[f][s].keys()
        #except:
        #    scaffold[s]=d[f][s].keys()
        for p in list(set(d[f][s].keys())):
            try:
                position[s][p]=d[f][s][p][0]
            except:
                position[s]={}
                position[s][p]=d[f][s][p][0]
seq=""
for s in position:
    for p in sorted(position[s]):
        seq+=position[s][p]
print ">ref" 
print seq
for f in d:
    print ">"+f
    seq=""
    #for s in scaffold:
    for s in position:
        for p in sorted(position[s]):
            try:
                seq+=d[f][s][p][1]
            except:
                seq+=position[s][p]
    print seq
