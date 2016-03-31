#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2010, 陈同 (chentong_biology@163.com).  
Please see the license file for legal information.
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
import sys
read=1
out=1

def main():
    print >>sys.stderr, "Print the result to screen"
    if len(sys.argv) < 3:
        print >>sys.stderr, 'Using python %s groupfile annofile [annofile]' \
            % sys.argv[0]
        sys.exit(0)
    
    #------------------------------------------------------
    global read
    global out
    annodict={}
    '''
    annodict={ locus : \n    \tdescription or
        interprodescription\nCC:\tGoid:\tGoterm\n
        MF:\tGoid:\tGoterm\n   \tGoid:\tGoterm\nBP:\tGoid:\tGoterm\n\n\n
             }
    '''
    list42 = (0, 3, 5, 1)
    dict42 = {3:"CC:", 5:"MF:", 1:"BP:"}
    #description0 1bp_id 2bp_name 3cc_id 4cc_name 5mf_id 6mf_name 7interpro 8interpro_description
    for file in sys.argv[2:]:
        header = 1
        locus = ''
        annolist = []
        for line in open(file):
            if header:
                header -= 1
            else:
                locus1, anno = line.split('\t',1)
                if locus1 != locus:
                    if locus:
                        if read:
                            print >>sys.stderr, locus
                        annodict[locus]='\n'
                        lenannolist = len(annolist)
                        for i42 in list42:
                            if not i42:
                                hasDes = 0
                                for i in range(lenannolist):
                                    if annolist[i][0]:
                                        annodict[locus] += annolist[i][0]+'\n'
                                        hasDes = 1
                                        break
                                interprodict = {}
                                for i in range(lenannolist):
                                    interproid = annolist[i][7]
                                    if interproid:
                                        hasDes = 1
                                        if interprodict.has_key(interproid):
                                            continue
                                        else:
                                            interprodict[interproid] =\
                                                annolist[i][8]
                                if len(interprodict):
                                    for keyvalue in interprodict.items():
                                        annodict[locus] += ':\t'.join(keyvalue)
                                        #annodict[locus] += ''.join(interprodict.values())
                                if not hasDes:
                                    annodict[locus] += 'No description\n'
                            #------------------------------------------------
                            else:
                                godict = {}
                                for i in range(lenannolist):
                                    gokey = annolist[i][i42]
                                    if gokey:
                                        if godict.has_key(gokey):
                                            continue
                                        else:
                                            godict[gokey]= '\t'+gokey+':\t'+annolist[i][i42+1]+'\n'
                                if len(godict):
                                    annodict[locus] += dict42[i42]+\
                                        ''.join(godict.values())
                            #---else------------------------------------------
                        #---for--i42 in list42--------------------------------
                    #-------if--locus-----------------------------------------
                    annolist = []
                    locus = locus1
                #-----------if--locus1 != locus----------------------------------
                annolist.append(anno.split('\t'))
            #-----else-----------------------------------------------------------
        #deal with the last one
        if read:
            print >>sys.stderr, locus
        annodict[locus]='\n'
        lenannolist = len(annolist)
        for i42 in list42:
            if not i42:
                hasDes = 0
                for i in range(lenannolist):
                    if annolist[i][0]:
                        annodict[locus] += annolist[i][0]+'\n'
                        hasDes = 1
                        break
                interprodict = {}
                for i in range(lenannolist):
                    interproid = annolist[i][7]
                    if interproid:
                        hasDes = 1
                        if interprodict.has_key(interproid):
                            continue
                        else:
                            interprodict[interproid] =\
                                annolist[i][8]
                if len(interprodict):
                    for keyvalue in interprodict.items():
                        annodict[locus] += ':\t'.join(keyvalue)
                    #annodict[locus] += ''.join(interprodict.values())
                if not hasDes:
                    annodict[locus] += 'No description\n'
            #------------------------------------------------
            else:
                godict = {}
                for i in range(lenannolist):
                    if annolist[i][i42]:
                        gokey = annolist[i][i42]
                        if gokey:
                            if godict.has_key(gokey):
                                continue
                            else:
                                godict[gokey]= '\t'+gokey+':\t'+annolist[i][i42+1]+'\n'
                if len(godict):
                    annodict[locus] += dict42[i42]+\
                        ''.join(godict.values())
            #---else------------------------------------------
        #---End with the last one-----------------------------
        #-------End reading one file----------------------------
        if 0:
            for key118, value118 in annodict.items():
                print key118
                print value118
    #-----------End reading multiple file--------------------
    #--------------------------------------------------------------------------------
    for line in open(sys.argv[1]):
        print line,
        if line[0].isalnum():
            locus = line.split('\t',2)[1]
            locus = locus.split()[0]
            if out:
                print >>sys.stderr, "-----------", locus
            if annodict.has_key(locus):
                print annodict[locus]
            else:
                print 'No this locusus'
                #print >>sys.stderr, 'No this locusus'
                #sys.exit(0)
    #--------------------------------------------------
if __name__ == '__main__':
    main()

