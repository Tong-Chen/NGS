#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2017, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:

    First column represents Molecular weight or amplicon size
    
    First row represents sample anmes

    Other values represenets concentrations in unit of ng/vl

    SIZE    samp1   samp2   samp3   samp4
    100   0   0   0   0
    200   1   0   0   0
    300   1   1   0   1
    400   0   0   1   0

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
import math
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0

def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
    print json_dumps(content,indent=1)

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="A matrix file with format specified above.")
    parser.add_option("-g", "--gel-concentration", dest="gel_concent",
        default=2, type="float", help="A number. Default <2> represents 2%.")
    parser.add_option("-v", "--voltage", dest="voltage",
        default=20, type="float", help="Voltage. Default <20> V.")
    parser.add_option("-t", "--time", dest="time",
        default=40, type="float", help="Run time. Default <40> minutes.")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def gamma(x, t, k=2):
    b = 1.0 / t
    x = x * k
    c = 0.0
    
    for i in range(0, k):
        c += (math.exp(-b * x) * (b * x) ** i) / math.factorial(i)

    return c
#-------END gamma-------------------
def readInGelFile(file):
    '''
    Molecular weight or amplicon size

    SIZE    samp1   samp2   samp3   samp4
    size1   0   0   0   0
    size2   200   0   0   0
    size3   200   200   0   200
    size4   0   0   200   0
    '''
    
    header = 1
    aDict = {}
    for line in open(file):
        if header:
            sampleL = line.split()
            for sample in sampleL[1:]:
                aDict[sample] = []
            header -= 1
            continue
        #_-----------------------
        lineL = line.split()
        key = float(lineL[0])
        for sample, concentration in zip(sampleL[1:], lineL[1:]):
            concentration = int(concentration)
            if concentration:
                aDict[sample].append([key,concentration])
        #----------------------------------------------
    return aDict,sampleL[1:]
#----------------------------


def plot(band_matrix, xpositionS, xvalueS, markerPositionS, markerSizeS):

    output = band_matrix + '.r'
    output_fh = open(output, 'w')
    print >>output_fh, """
usePackage <- function(p) {{
	    if (!is.element(p, installed.packages()[,1]))
			        install.packages(p, dep = TRUE)
 	   require(p, character.only = TRUE)
}}

usePackage("ggplot2")
usePackage("reshape2")


data <- read.table(file="{file}", sep="\t", header=T, row.names=1,
	check.names=F, quote="")


data$id <- rownames(data)
idlevel <- as.vector(rownames(data))

idlevel <- rev(idlevel)

data.m <- melt(data, c("id"))

data.m$id <- factor(data.m$id, levels=idlevel, ordered=T)

p <- ggplot(data=data.m, aes(x=variable, y=id)) + 	geom_tile(aes(fill=value)) 

midpoint = 0

p <- p + scale_fill_gradient2(low="black", mid="grey",
	high="white", midpoint=midpoint, name=" ",
	na.value="grey")

p <- p + theme(axis.ticks=element_blank()) + theme_bw() + 
	theme(panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), panel.border = element_blank()) + xlab("") + 
	ylab("") + labs(title="")

p <- p + theme(axis.text.x=element_text(angle=45,hjust=0, vjust=0))

none='none'
legend_pos_par <- none

p <- p + theme(legend.position=legend_pos_par)


xtics_pos <- {xpositionS}
xtics_value <- {xvalueS}

p <- p + scale_x_discrete(breaks=xtics_pos, labels=xtics_value, position="top")

ytics_pos <- {markerPositionS}
ytics_value <- {markerSizeS}

p <- p + scale_y_discrete(breaks=ytics_pos, labels=ytics_value, position="right")

p <- p + theme(text=element_text(size=14))


ggsave(p, filename="{file}.png", dpi=300, width=10, height=20, units=c("cm"))
""".format(file=band_matrix, xvalueS=xvalueS, xpositionS=xpositionS, markerSizeS=markerSizeS, 
        markerPositionS=markerPositionS)

    output_fh.close()
    os.system("Rscript "+output)

#----------------------plot----------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    #verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    
    gel_concentration = options.gel_concent
    optimum_DNA_length = 2000 / gel_concentration ** 3

    # samples: {sampleName: [(expected DNA length, DNA concentration), (expected DNA length2, DNA concentration2)],)
    # sampleD = {"A": [(400,200)],"B":[(500,200)],"C":[(900,200), (1000,200)]}
    sampleD, nameL = readInGelFile(file)

    # Add ladder
    nameL.append("Markers")
    sampleD["Markers"] = [(50,200), (100,200),(150,200),(200,200),(250,200),(300,200),(350,200),(400,200),(500,200),(1000,200)]

    len_sampleD = len(sampleD)


    # Number of lanes
    lane_count = len_sampleD
    # Width and height of dingle lane
    lane_width = 30
    lane_height = 3
    # Intervals between neighboring lanes
    lane_interval = 6
    gel_border  = 12

    gel_width = 2 * gel_border + lane_count * lane_width + (lane_count-1) * lane_interval

    #gel_height = gel_width

    # Generate empty lanes
    # Every 1 unit from up_border to all height
    strandD = {}


    # Loading samples
    # position: represents start position, loaing hole
    for name in nameL:
        strandD[name] = []
        dnaL = sampleD[name]
        for dna in dnaL:
            band = {"size": float(dna[0]), 'conc': dna[1]*1.0/lane_height, 'position': gel_border+5}
            strandD[name].append(band)
    
    # Run 
    # Move loaded DNA down the gel at a rate dependent on voltage, DNA length and agarose concentration.
    time = options.time
    voltage = options.voltage
    max_dist = 0.25 * time * voltage

    maxposition = 0
    for name in nameL:
        for bandD in strandD[name]:
            g = gamma(bandD['size']/20, int(optimum_DNA_length/20))   
            bandD['position'] += max_dist * g
            if maxposition < bandD['position']:
                maxposition = bandD['position']
            #print >>sys.stderr, bandD
    #-----------------------------------------------------
    maxposition = int(maxposition+gel_border+1)
    # Determines where in the concentration of DNA in every part of the gel
    # Generate an empty lane, O represents no band

    markerSize     = [str(int(bandD["size"])) for bandD in strandD["Markers"]]
    markerSize.reverse()
    markerSizeS = "c("+','.join(["'"+i+"'" for i in markerSize])+ ")"
    markerPosition = ['CT'+str(int(bandD["position"]+0.5)) for bandD in strandD["Markers"]]
    markerPosition.reverse()
    markerPositionS = "c("+','.join(["'"+i+"'" for i in markerPosition])+ ")"
    #print >>sys.stderr, markerSizeS
    #print >>sys.stderr, markerPositionS
    #laneL = [[-1 for j in range(maxposition)] for i in range(len_sampleD)]

    laneL = []
    for i in range(len_sampleD):
        tmpL = [-1000 for j in range(maxposition)]
        tmpL[gel_border] = -10
        tmpL[gel_border-1] = 0
        tmpL[gel_border-2] = -10
        laneL.append(tmpL)
    #-------------------------------------
    band_count = maxposition

    for i in range(len_sampleD):
        name = nameL[i]
        for bandD in strandD[name]:
            for y in range(lane_height-2):
                pos = int(bandD['position'])+y
                if pos < band_count - 4:
                    laneL[i][pos-1] += 0.12 * bandD['conc'] * bandD['size']
                    laneL[i][pos]   += 0.2  * bandD['conc'] * bandD['size']
                    laneL[i][pos+1] += 0.12 * bandD['conc'] * bandD['size']
                #-----Blur edges-----------------------------
            #---------------------------
        #_--------------------------
    #-------------------------------
    max_value = max([max(i) for i in laneL])
    min_value = -1 * max_value
    for i in range(len_sampleD):
        for j in range(maxposition):
            if laneL[i][j] == -1000:
                laneL[i][j] = min_value
    #print laneL    
    # Expose
    # print "ID\t{}".format("\t".join(nameL))
    newCol = len(nameL) * 8 + 2
    xposition = ['ct'+str(i) for i in range(3, newCol,8)]

    xpositionS = "c("+','.join(["'"+i+"'" for i in xposition])+ ")"
    xvalueS = "c("+','.join(["'"+i+"'" for i in nameL])+ ")"
    #print >>sys.stderr, xpositionS
    #print >>sys.stderr, '\t'.join(nameL)

    # Plot data

    band_matrix = file + '.virtualGel.data'
    band_matrix_fh = open(band_matrix, 'w')
    print >>band_matrix_fh, "ID\t{}".format('\t'.join(['ct'+str(i) for i in range(newCol)]))
    for i in range(maxposition):
        tmpL = ["CT"+str(i)]
        tmpL.append(str(min_value))
        tmpL.append(str(min_value))
        for lane in laneL:
            # Blur edges, 6 main lane, two margin lane
            if lane[i] > 0:
                tmpL.append(str(0.8*lane[i]))
                tmpL.append(str(0.9*lane[i]))
            else:
                tmpL.append(str(lane[i]))
                tmpL.append(str(lane[i]))
            tmpL.append(str(lane[i]))
            tmpL.append(str(lane[i]))

            if lane[i] > 0:
                tmpL.append(str(0.9*lane[i]))
                tmpL.append(str(0.8*lane[i]))
            else:
                tmpL.append(str(lane[i]))
                tmpL.append(str(lane[i]))

            tmpL.append(str(min_value))
            tmpL.append(str(min_value))
        print >>band_matrix_fh, "\t".join(tmpL)
    band_matrix_fh.close()
    # Generate R script for plot 
    
    plot(band_matrix, xpositionS, xvalueS, markerPositionS, markerSizeS)

    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    #if verbose:
    #    print >>sys.stderr,            "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " %         (' '.join(sys.argv), startTime, endTime)
    fh.close()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------
#

