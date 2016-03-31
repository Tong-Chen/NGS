#!/usr/bin/python
import sys
if len(sys.argv) < 2:
    print >>sys.stderr, "%s Fasta_file DNA/RNA ignor_nonATCG[default no transfer. If \
0 given, it will delete other characters]" % sys.argv[0]
    sys.exit(1)
fasta = sys.argv[1]
if len(sys.argv) > 2:
    type = sys.argv[2]
else:
    type = 'DNA'

if len(sys.argv) >3:
    ignore = sys.argv[3]
else:
    ignore = 1

#--------------------------------------
if type == 'DNA':
    dict = {'A':'T', 'G':'C', 'T':'A','C':'G','a':'t','g':'c','t':'a','c':'g'}
elif type == 'RNA':
    dict = {'A':'U', 'G':'C', 'U':'A','C':'G','a':'u','g':'c','u':'a','c':'g'}
#-------------------------------------------

seqL = []

for line in open(fasta):
    if line[0] == '>':
        if len(seqL) > 0:
            seq = ''.join(seqL)
            seqL = []
            strarray = list(seq)
            strarray.reverse()
            reverseCom = []

            for i in strarray:
                if i not in dict:
                    if ignore:
                        reverseCom.append(i)
                else:
                    reverseCom.append(dict[i])
            #-----------------------------------
            print key
            print ''.join(reverseCom)
        #---------------------------------------
        key = line.strip()
    else:
        seqL.append(line.strip())
#------------------------------------------------
if len(seqL) > 0:
    seq = ''.join(seqL)
    seqL = []
    strarray = list(seq)
    strarray.reverse()
    reverseCom = []

    for i in strarray:
        if i not in dict:
            if ignore:
                reverseCom.append(i)
        else:
            reverseCom.append(dict[i])
    #-----------------------------------
    print key
    print ''.join(reverseCom)

#---------------------------------------
