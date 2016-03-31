#!/usr/bin/env python

#Last Change:  2010-11-03 13:38:36
def compute_prefix_function(p):
    '''
    compute_prefix_function(p) --> a list

    [p] is a pattern sequence to be searched in a longer sequence.
    '''
    m = len(p)
    pi = [0] * m
    k = 0
    for q in range(1, m):
        while k > 0 and p[k] != p[q]:
            k = pi[k - 1]
        if p[k] == p[q]:
            k = k + 1
        pi[q] = k
    return pi

def kmp_matcher(t, p):
    '''
    kmp_matcher(t, p) --> a list

    [t, p] is two sequences, which 'p' maybe partof 't'.
    '''
    n = len(t)
    m = len(p)
    pi = compute_prefix_function(p)
    q = 0
    index = []
    for i in range(n):
        while q > 0 and p[q] != t[i]:
            q = pi[q - 1]
        if p[q] == t[i]:
            q = q + 1
        if q == m:
            q = 0
            index.append(i-m+2)
    #-----end of for---------------------------
    return -1 if len(index) == 0 else index

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print >>sys.stderr, "%s string substr" % sys.argv[0]
        sys.exit(1)
    str = sys.argv[1]
    substr = sys.argv[2]
    #str = "1230012123211212213012300121232112122130123001212321121221312300121232112122130"
    #substr = "1212"
    #for i in range(10000000):
    #    kmp_matcher(str, substr)
    print kmp_matcher(str, substr)
    #print compute_prefix_function('ababacb')
