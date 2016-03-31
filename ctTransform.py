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

codonTable = {
    'TCA' : 'S',    # Serine 
    'TCC' : 'S',    # Serine 
    'TCG' : 'S',    # Serine 
    'TCT' : 'S',    # Serine 
    'TTC' : 'F',    # Phenylalanine 
    'TTT' : 'F',    # Phenylalanine 
    'TTA' : 'L',    # Leucine 
    'TTG' : 'L',    # Leucine 
    'TAC' : 'Y',    # Tyrosine 
    'TAT' : 'Y',    # Tyrosine 
    'TAA' : '_',    # Stop 
    'TAG' : '_',    # Stop 
    'TGC' : 'C',    # Cysteine 
    'TGT' : 'C',    # Cysteine 
    'TGA' : '_',    # Stop 
    'TGG' : 'W',    # Tryptophan 
    'CTA' : 'L',    # Leucine 
    'CTC' : 'L',    # Leucine 
    'CTG' : 'L',    # Leucine
    'CTT' : 'L',    # Leucine 
    'CCA' : 'P',    # Proline 
    'CCC' : 'P',    # Proline 
    'CCG' : 'P',    # Proline 
    'CCT' : 'P',    # Proline 
    'CAC' : 'H',    # Histidine 
    'CAT' : 'H',    # Histidine 
    'CAA' : 'Q',    # Glutamine 
    'CAG' : 'Q',    # Glutamine 
    'CGA' : 'R',    # Arginine 
    'CGC' : 'R',    # Arginine 
    'CGG' : 'R',    # Arginine 
    'CGT' : 'R',    # Arginine 
    'ATA' : 'I',    # Isoleucine 
    'ATC' : 'I',    # Isoleucine 
    'ATT' : 'I',    # Isoleucine 
    'ATG' : 'M',    # Methionine 
    'ACA' : 'T',    # Threonine 
    'ACC' : 'T',    # Threonine 
    'ACG' : 'T',    # Threonine 
    'ACT' : 'T',    # Threonine 
    'AAC' : 'N',    # Asparagine 
    'AAT' : 'N',    # Asparagine 
    'AAA' : 'K',    # Lysine 
    'AAG' : 'K',    # Lysine 
    'AGC' : 'S',    # Serine 
    'AGT' : 'S',    # Serine 
    'AGA' : 'R',    # Arginine 
    'AGG' : 'R',    # Arginine 
    'GTA' : 'V',    # Valine 
    'GTC' : 'V',    # Valine 
    'GTG' : 'V',    # Valine 
    'GTT' : 'V',    # Valine 
    'GCA' : 'A',    # Alanine 
    'GCC' : 'A',    # Alanine 
    'GCG' : 'A',    # Alanine 
    'GCT' : 'A',    # Alanine 
    'GAC' : 'D',    # Aspartic Acid 
    'GAT' : 'D',    # Aspartic Acid 
    'GAA' : 'E',    # Glutamic Acid 
    'GAG' : 'E',    # Glutamic Acid 
    'GGA' : 'G',    # Glycine 
    'GGC' : 'G',    # Glycine 
    'GGG' : 'G',    # Glycine 
    'GGT' : 'G',    # Glycine 
    }

def translateEasy(seq, dna=1):
    if not dna:
        seq = seq.replace('U','T')
    pro = ''
    lenseq = len(seq)
    for i in range(0,lenseq,3):
        pro += codonTable[seq[i:i+3]]
    return pro
#----------------------------------------

if __name__ == '__main__':
    print >>sys.stderr, 'A module'
