#!/usr/bin/env python

#script to identify transcripts with an extended 3'UTR post-MBT relative to pre-MBT
#Input is a bed file with the last exons of all loci (where multiple last exons; the exon furthest away from the TSS)
#and bamfiles from pre and post-MBT

import sys, subprocess

bed = open(sys.argv[1],'r').readlines()
bamfile_pre = sys.argv[2] 
bamfile_post = sys.argv[3]

for i in bed:
	tmp = i.split()
	chr = tmp[0].replace('chr','')
	start = int(tmp[1])
	end = int(tmp[2])
	length = int(end) - int(start)
	half = int(round(length/2.0))
	parta_s = int(start)
	parta_e = int(start + half)
	partb_s =int(parta_e + 1)
	partb_e = int(partb_s + half)
	p1 = str(chr+':'+str(parta_s)+'-'+str(parta_e))
	p2 = str(chr+':'+str(partb_s)+'-'+str(partb_e))		
	exon = 'chr'+str(chr)+':'+str(start)+'-'+str(end)
	strand = tmp[5]
	#print p1
	#print p2

	#Run samtools and wc -l to get the number of reads
	pre_p1 = subprocess.Popen(['samtools','view',bamfile_pre, p1], stdout = subprocess.PIPE)	
	pre_p2 = subprocess.Popen(['samtools','view',bamfile_pre, p2], stdout = subprocess.PIPE)
	pre_p1_c = subprocess.Popen(['wc','-l'], stdin = pre_p1.stdout, stdout = subprocess.PIPE)
	pre_p2_c = subprocess.Popen(['wc','-l'], stdin = pre_p2.stdout, stdout = subprocess.PIPE)	
	pre_p1_co = int(pre_p1_c.communicate()[0])
	pre_p2_co = int(pre_p2_c.communicate()[0])

	
	post_p1 = subprocess.Popen(['samtools','view',bamfile_post, p1], stdout = subprocess.PIPE)
        post_p2 = subprocess.Popen(['samtools','view',bamfile_post, p2], stdout = subprocess.PIPE)
	post_p1_c = subprocess.Popen(['wc','-l'], stdin = post_p1.stdout, stdout = subprocess.PIPE)
        post_p2_c = subprocess.Popen(['wc','-l'], stdin = post_p2.stdout, stdout = subprocess.PIPE)
        post_p1_co = int(post_p1_c.communicate()[0])
        post_p2_co = int(post_p2_c.communicate()[0])

	#calculate score
	#print pre_p1_co,pre_p2_co,post_p1_co,post_p2_co

	try:
		q1 = (pre_p2_co+1)/(float(pre_p2_co)+float(pre_p1_co)+1)
		q2 = (post_p2_co+1)/(float(post_p2_co+1)+float(post_p1_co)+1)
		out = q1/q2
		tot_pre = float(pre_p2_co)+float(pre_p1_co)
		tot_post = float(post_p2_co+1)+float(post_p1_co)

		print tmp[3], exon,strand,p1,p2,pre_p1_co,pre_p2_co,post_p1_co,post_p2_co,out

	except:
		print tmp[3],exon, 'Problem?'
		

