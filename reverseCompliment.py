#!/usr/bin/python
import sys
if len(sys.argv) < 3:
    print >>sys.stderr, "%s seq type[DNA/RNA] ignor_nonATCG[default no transfer \
(a value represents true). If 0 given, it will delete other characters]" % sys.argv[0]
    sys.exit(1)
str = sys.argv[1]
type = sys.argv[2]
if len(sys.argv) > 3:
    ignore = sys.argv[3]
else:
    ignore = 1
#--------------------------------------
if type == 'DNA':
    dict = {'A':'T', 'G':'C', 'T':'A','C':'G','a':'t','g':'c','t':'a','c':'g'}
elif type == 'RNA':
    dict = {'A':'U', 'G':'C', 'U':'A','C':'G','a':'u','g':'c','u':'a','c':'g'}


strarray = list(str)
strarray.reverse()
reverseCom = []

for i in strarray:
    if i not in dict:
        if ignore:
            reverseCom.append(i)
    else:
        reverseCom.append(dict[i])

print ''.join(reverseCom)
