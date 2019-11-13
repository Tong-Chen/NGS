#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2017, 陈同 (chentong_biology@163.com),柴国师 (cgs@immunet.cn) 
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This script is designed to generate a circos figure by supplying the chrmosomes bed file, bands bed file(optional), a or a few highlights bed files, a or a few heatmap bed files, a name of the circos.conf file. Following is the data format of these files.

1. **chrmosomes bed file**

chromosomes.sizes.bed ('#' started line will be skipped, only used for description):

## <main ID column>: 
    # For single specie circos plot, the first column will be used as main ID column.
    # For multiple species circos plot, the forth column will be used as main ID column. 
    # Main ID column represents the names of current regions. All sub-regions must have their parent fields same as ID columns.
## 
## <label column> (the forth column) will be the label text along circles
## 

#chromosomeID	start	end	label
chr1	0	248956422	hs1
chr2	0	242193529	hs2
chr3	0	198295559	hs3
chr4	0	190214555	hs4
chr5	0	181538259	hs5

2. bands bed file (optional)

GRCh38.cytoBand.bed ('#' started line will be skipped, only used for description)

## <ParentID>: 
    # For multiple species plot, the last column will be treated as Parent ID column, which should be same as label column in <chrmosomes bed file>. 
    # For single specis plot, first column will be treated as Parent ID.

#chromosomeID	start	end	arm	color   ParentID
chr1	0	2300000	p36.33	gneg    hs1
chr1	2300000	5300000	p36.32	gpos25  hs1
chr1	5300000	7100000	p36.31	gneg    hs1
chr1	7100000	9100000 p36.23	gpos25  hs1
chr1	9100000	12500000	p36.22	gneg    hs1


3. **highlights bed file**

GRCh38.genome_landmark.TSS_UP1k.bed('#' started line will be skipped,  only used for description):


## <ParentID>: 
    # For multiple species plot, the last column will be treated as Parent ID column, which should be same as label column in <chrmosomes bed file>. 
    # For single specis plot, first column will be treated as Parent ID.
## Other columns are standard bed files. 
## The forth column would be used as labels to show along plots.

#chromosomeID	start	end	label   ParentID
chr1	181392	182392	TSS_UP1k    hs1
chr1	200321	201321	TSS_UP1k    hs1
chr1	922927	923927	TSS_UP1k	hs1
chr1	959308	960308	TSS_UP1k	hs1
chr1	959586	960586	TSS_UP1k	hs1
chr1	965496	966496	TSS_UP1k	hs1

4. **heatmap bed file**

Af_Z101_C.CpG.met.R8.1M.bed ('#' started line will be skipped,  only used for description):

## <ParentID>: 
    # For multiple species plot, the last column will be treated as Parent ID column, which should be same as label column in <chrmosomes bed file>. 
    # For single specis plot, first column will be treated as Parent ID.
## Other columns are standard bed files. 
## The forth column would be used as values to display in heatmap

chromosomeID	start	end	values  ParentID
chr4	0	1000000	4310.19999999999708962  hs4
chr4	1000000	2000000	6920.70000000000163709  hs4
chr4	2000000	3000000	3210.05999999999085048  hs4
chr4	3000000	4000000	4189.57000000002517481  hs4
chr4	4000000	5000000	940.520000000001573426  hs4
chr4	5000000	6000000	768.859999999998990461  hs4
chr4	6000000	7000000	2072.13000000000511136  hs4

The name of the circos.conf file can be circos1.conf, circos2.conf and so on. The name must be ended with .conf.

cmd:
circos.py -c chromesomes.bed -b GRCh38.cytoBand.bed -H "GRCh38.genome_landmark.TSS_UP1k_1.bed GRCh38.genome_landmark.TSS_UP1k_2.bed GRCh38.genome_landmark.TSS_UP1k_3.bed" -m "Af_Z101_C.CpG.met.R8.1M.bed Au_113_F.CpG.met.R8.1M.bed Nu_E123_D.CpG.met.R8.1M.bed" -n cgs1 -s 2.0


'''

import sys
import re
import os
import random
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-M", "--multiple-species", dest="multi_spe", 
        default=False, action="store_true", help="Specify whether single specie or multiple species are used for plotting. For multi_spe, input format is a little complex.")
    parser.add_option("-c", "--chrmosomes", dest="chrmosomes", help="Supply a chrmosomes bed file with format specified above.")
    parser.add_option("-b", "--bands", dest="bands", help="This file is optional, a bed file with format specified above. Currently this will be ignored. One can only set bands using <-k>")
    parser.add_option("-k", "--known-species", dest="known_spe", help="Specify the genome version to get specific karyotype files saved in CIRCOS. Currently support <hg19, hg38, tair10, pt4, dm3, mm10, mm9, rn4, yeast, zeamays, sorghum, oryzasativa>. This parameter will make <-c> and <-b> ignored.")
    parser.add_option("-H", "--highlights", dest="highlights", help="Supply one or a few bed files. The first field is the same as the field in the chrmosomes bed file and band bed file. Second and third fields are the start and end of functional elements such as transcription start site(TSS). The forth field is the name of functional elements. The last column would be Parent ID when <-m> specified,  other wise the first column would be Parent ID. Detailed format above.")
    parser.add_option("-v", "--valueFile", dest="heatmap", help="Supply one or a few bed files. The first field is the same as the field in the chrmosomes bed file and band bed file and highlights bed file. Second and Third fields are the start and end of a value such as the level of Methylation site. The fourth field is the value of Methylation site. The last column would be Parent ID when <-m> specified, otherwise the first column would be Parent ID. Detailed format above.")
    parser.add_option("-t", "--type", default='heatmap', dest="plot_type", help="Plot type for given value files. Default <heatmap>. Accept <line>,  <histogram>, <scatter>. If only one value is given, all value file will be plotted in same format. If multiple values are given, each will have specified plot types.")
    parser.add_option("-n", "--filename", dest="circos", help="Supply a prefix of the circos.conf file name and the circos figure name. The files will be created by this program.")
    parser.add_option("-s", "--scale", dest="scale", default="0.5", help="The value of scale_log_base. This argument is optional, and the default value is 0.5. The values of scale_log_base can be 0.5, 1.0, 2.0. 0.5 means greater dynamic range of color for smaller values; 1.0 means colors uniformly distributed across range of values; 2.0 means greater dynamic range of color for larger values.")
    (options, args) = parser.parse_args(argv[1:])
    assert options.chrmosomes != None or options.known_spe != None, "A filename needed for -c"
    assert options.circos != None, "A filename needed for -n"
    return (options, args)
#--------------------------------------------------------------------

# 生成circos支持的染色体文件
def outfile_chromosomes(chrom, bands, chromosomes_file, multi_spe):
    chromosomes_color_list=["chr1","chr2","chr3","chr4","chr5","chr6",
        "chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
        "chr15","chr16","chr17","chr18","chr19","chr20","chr21",
        "chr22","chrx","chry"]
    chromosomes_color_dict = dict([(i.upper(), i) for i in chromosomes_color_list])

    # save colors
    colorD = {}
    out = open(chromosomes_file,"w")
    parentIDl = []
    for line in open(chrom):
        if line[0] == '#':
            continue
        lineL = line.strip().split("\t")
        if multi_spe:
            parentID = lineL[3]
        else:
            parentID = lineL[0]
        parentIDl.append(parentID)
        chr_name = lineL[0].upper()
        color = colorD.get(chr_name, 
            chromosomes_color_dict.get(
                chr_name, 
                random.sample(chromosomes_color_list, 1)[0])
            )
        colorD[chr_name] = color

        tmpL = ['chr', '-', parentID, lineL[3], lineL[1], lineL[2], color]
        print >>out, ' '.join(tmpL)

    #if bands:
    #    for line in open(bands):
    #        line = line.strip()
    #        list_tmp1=i.split("\t")
    #        if list_tmp1[0] in chrmosomes_list:
    #            out.write("band "+list_tmp1[0]+" "+list_tmp1[3]+" "+list_tmp1[3]+" "+list_tmp1[1]+" "+list_tmp1[2]+" "+list_tmp1[4]+"\n")
    out.close()
    return parentIDl
#------------------------------------------------------



def outfile_bands(file1,file4):
    out=open(file1+".txt","a")
    fh=open(file4)
    for i in fh:
        i=i.strip()
        list_tmp1=i.split("\t")
        out.write("band "+list_tmp1[0]+" "+list_tmp1[3]+" "+list_tmp1[3]+" "+list_tmp1[1]+" "+list_tmp1[2]+" "+list_tmp1[4]+"\n")
    out.close()

# 生成circos支持的高亮元件的文件
def outfile_highlights(highlights_list, multi_spe, circos):
    if multi_spe:
        id_index = -1
    else:
        id_index = 0
    #--------------------------------
    colorD = {}

    highlights_color_list=["red","green","blue","purple","yellow",
        "orange","grey", "lpred", "lpblue", "lpurple", "lporange", 
        "dpurple", "dporange", "dyellow", "dpgreen", "dpred" 
        "vdpred", "vdpblue","vdpurple", "vdporange", "vdyellow",
        "vdgreen", "vvlpred", "vvlpblue", "vvlppurple", "vvlporange", 
        "vvlyellow", "vvlpgreen", "vlpblue", "vlpgreen", "vlporange", 
        "vlppurple", "vlpred", "vlyellow"]

    highlightL = []
    for file in highlights_list:
        output = file + '.circos_input.txt'
        highlightL.append(output)
        output_fh = open(output, 'w')
        for line in open(file):
            if line[0] == '#':
                continue
            lineL = line.strip().split('\t')
            ID = lineL[id_index]
            start = lineL[1]
            end = lineL[2]
            label = lineL[3]
            #fill = colorD.get(label, highlights_color_list.pop())
            fill = colorD.get(label, 'unknown')
            if fill == 'unknown':
                fill = highlights_color_list.pop()
                colorD[label] = fill
            tmpL = [ID, start, end, "fill_color="+fill]
            print >>output_fh, ' '.join(tmpL)
        output_fh.close()
    #----------------------------------------------------------------
    #high_legend = circos + '.highlight.legend.txt'
    #high_legend_fh = open(high_legend, 'w')
    #for label, color in colorD:
    #    print >>high_legend_fh, '\t'.join([label, color])
    #high_legend_fh.close()

    return highlightL, colorD
#------------------------------------------------------------------------

# 生成circos支持的热图数据文件

def outfile_heatmap(inputHeatmapL, multi_spe):
    if multi_spe:
        id_index = -1
    else:
        id_index = 0
    #--------------------------------
    heatmapL = []
    min_v = 100000000
    max_v = -100000
    for file in inputHeatmapL:
        output = file + '.circos_input.txt'
        heatmapL.append(output)
        output_fh = open(output, 'w')
        for line in open(file):
            lineL = line.strip().split('\t')
            ID = lineL[id_index]
            lineL[0] = ID
            value = float(lineL[3])
            if min_v > value:
                min_v = value
            if max_v < value:
                max_v = value
            print >>output_fh, ' '.join(lineL[:4])
        output_fh.close()
    #----------------------------
    return heatmapL, min_v, max_v
#-------------------------------------------------------


# 生成circos.conf文件
def circos_header(out):
    print >>out, "<<include etc/colors_fonts_patterns.conf>>"
    print >>out, "<<include ideogram.conf>>"
    print >>out, "<<include ticks.conf>>\n"

def circos_image(out,image_name,image_position="./"):
    out.write("<image>\n<<include etc/image.conf>>\n")
    out.write("file*="+image_name+"\n")
    out.write("dir*="+image_position+"\n</image>\n")

def circos_chromosomes_display(out,chromosome_file):
    out.write("karyotype = "+chromosome_file+"\n")
    out.write("chromosomes_units =1000000\n")
    out.write("chromosomes_display_default = yes\n")


def heatmap_scale(out,scale="0.5"):
    #out=open(fileX,"a")
    out.write("scale_log_base="+scale+"\n")
    #out.close()

#def circos_heatmap_display(out,file_num):
#    #out=open(fileX,"a")
#    list1=[]
#    dict2=dict()
#    span_num=0.6/file_num
#    for i in range(file_num):
#        list1.append([0.4+i*span_num,0.4+(i+1)*span_num])
#    list1[-1]=1.0
#    for j in range(len(list1)):
#        dict2[list2[j]]=list1[j]
#    for key in dict2:
#        out.write("<plot>\n")
#        out.write("file=./"+key+".txt\n")
#        out.write("r0="+dict2[key][0]+"\nr1="+dict2[key][1]+"\n</plot>\n")
#    #out.close()

def ideogram_conf(path1="."):
    out=open(path1+"/ideogram.conf","w")
    out.write("<ideogram>\n<spacing>\ndefault = 0.005r\n</spacing>\n<<include ideogram.position.conf>>\n<<include ideogram.label.conf>>\n<<include bands.conf>>\nradius* = 0.90r\n</ideogram>\n")
    out.close()

def ideogram_position_conf(path1="."):
    out=open(path1+"/ideogram.position.conf","w")
    out.write("radius = 0.775r\nthickness = 30p\nfill = yes\nfill_color = black\nstroke_thickness = 2\nstroke_color = black\n")
    out.close()

def ideogram_label_conf(path1="."):
    out=open(path1+"/ideogram.label.conf","w")
    out.write("show_label = yes\nlabel_font = default\nlabel_radius = dims(ideogram,radius_outer) + 0.05r\nlabel_size = 24\nlabel_parallel   = yes\nlabel_case = upper\n")
    out.close()

def bands_conf(path1="."):
    out=open(path1+"/bands.conf","w")
    out.write("show_bands = yes\nfill_bands = yes\nband_stroke_thickness = 2\nband_stroke_color = white\nband_transparency = 0\n")
    out.close()

def ticks_conf(path1="."):
    out=open(path1+"/ticks.conf","w")
    out.write("show_ticks = yes\nshow_tick_labels = yes\n<ticks>\nradius = dims(ideogram,radius_outer)\norientation = out\nlabel_multiplier = 1e-6\ncolor = black\nsize = 20p\nthickness = 3p\nlabel_offset = 5p\nformat = %d\n<tick>\nspacing = 1u\nshow_label = no\nsize  = 10p\n</tick>\n<tick>\nspacing = 5u\nshow_label = no\nsize  = 15p\n</tick>\n<tick>\nspacing = 10u\nshow_label = yes\nlabel_size = 10p\n</tick>\n</ticks>\n")
    out.close()

def main():
    options, args = cmdparameter(sys.argv)
    multi_spe = options.multi_spe
    chrmosomes = options.chrmosomes
    bands = options.bands
    known_spe = options.known_spe
    highlights=options.highlights
    if highlights:
        highlights_list=highlights.replace(',', ' ').split()
    else:
        highlights_list=[]
    heatmap=options.heatmap
    plot_type = options.plot_type
    if heatmap:
        heatmap_list=heatmap.replace(',', ' ').split()
        plot_typeL = plot_type.replace(',', ' ').split()
        if len(plot_typeL) == 1:
            plot_typeL = plot_typeL * len(heatmap_list)
    else:
        heatmap_list=[]
    circos=options.circos
    scale_log_base = options.scale
    heatmap_highlights_list=highlights_list+heatmap_list

    chromosomes_file = circos + '.chromsomes.circos_input.txt'
    known_speD = {
        "tair10": "data/karyotype/karyotype.arabidopsis.tair10.txt",
        "arabidopsis": "data/karyotype/karyotype.arabidopsis.txt",
        "pt4": "data/karyotype/karyotype.chimp.pt4.txt",
        "chimp": "data/karyotype/karyotype.chimp.txt",
        "dm3": "data/karyotype/karyotype.drosophila.hires.dm3.txt",
        "dm3": "data/karyotype/karyotype.drosophila.lowres.dm3.txt",
        "drosophila": "data/karyotype/karyotype.drosophila.txt",
        "hg16": "data/karyotype/karyotype.human.hg16.txt",
        "hg17": "data/karyotype/karyotype.human.hg17.txt",
        "hg18": "data/karyotype/karyotype.human.hg18.txt",
        "hg19": "data/karyotype/karyotype.human.hg19.txt",
        "hg38": "data/karyotype/karyotype.human.hg38.txt",
        "human": "data/karyotype/karyotype.human.txt",
        "mm10": "data/karyotype/karyotype.mouse.mm10.txt",
        "mm9": "data/karyotype/karyotype.mouse.mm9.txt",
        "mouse": "data/karyotype/karyotype.mouse.txt",
        "oryzasativa": "data/karyotype/karyotype.oryzasativa.txt",
        "rn4": "data/karyotype/karyotype.rat.rn4.txt",
        "rat": "data/karyotype/karyotype.rat.txt",
        "3": "data/karyotype/karyotype.rm.3.txt",
        "rm": "data/karyotype/karyotype.rm.txt",
        "sorghum": "data/karyotype/karyotype.sorghum.txt",
        "yeast": "data/karyotype/karyotype.yeast.txt",
        "zeamays": "data/karyotype/karyotype.zeamays.txt"
    }
    if known_spe:
        chromosomes_file = known_speD[known_spe]
    else:
        chrmosomes_list = outfile_chromosomes(chrmosomes, bands, 
                chromosomes_file, multi_spe)
    

    if highlights_list:
        highlightL, colorD = outfile_highlights(highlights_list, multi_spe, circos)
    else:
        highlightL = []
        colorD = {}

    if heatmap_list:
        heatmapL, minv, maxv = outfile_heatmap(heatmap_list, multi_spe)
    else:
        heatmapL = []
        minv = maxv = 'ct'

    circos_conf = circos+".circos.conf"
    circos_conf_fh = open(circos_conf, 'w')
    circos_header(circos_conf_fh)
    circos_image(circos_conf_fh,circos+".circos.png")
    circos_chromosomes_display(circos_conf_fh, chromosomes_file)

    # 从0.4开始放置heatmap，highlight
    start = 0.4
    if heatmap_highlights_list:
        file_num = len(heatmap_highlights_list)
        span_num = (1-start) / file_num
        positionsL = []
        for i in range(file_num):
            positionsL.append([start+i*span_num,start+(i+0.9)*span_num])
        positionsL.reverse()

        count=0
        if highlightL:
            print >>circos_conf_fh, "<highlights>\n\n"
            for i in highlightL:
                circos_conf_fh.write("<highlight>\n")
                circos_conf_fh.write("file="+i+"\n")
                circos_conf_fh.write("r0="+str(positionsL[count][0])+"r\n")
                circos_conf_fh.write("r1="+str(positionsL[count][1])+"r\n")
                circos_conf_fh.write("</highlight>\n\n")
                count=count+1
            circos_conf_fh.write("</highlights>\n\n")
        if heatmapL:
            circos_conf_fh.write("<plots>\n\ncolor = spectral-7-div-rev\nstroke_thickness = 1\nstroke_color = black\nscale_log_base = "+scale_log_base+"\n")
            for i, j in zip(heatmapL, plot_typeL):
                circos_conf_fh.write("\n<plot>\n")
                circos_conf_fh.write("file="+i+"\n")
                circos_conf_fh.write("type="+j+"\n")
                circos_conf_fh.write("r0="+str(positionsL[count][0])+"r\n")
                circos_conf_fh.write("r1="+str(positionsL[count][1])+"r\n")
                circos_conf_fh.write("</plot>\n\n")
                count=count+1
            circos_conf_fh.write("</plots>\n")
        #-----------------------------
    circos_conf_fh.write("<<include etc/housekeeping.conf>>\nmax_points_per_track*  =  2500000\ndata_out_of_range* = trim\n")
    circos_conf_fh.close()
    ideogram_conf()
    ideogram_position_conf()
    ideogram_label_conf()
    bands_conf()
    ticks_conf()
    cmd="circos -conf "+circos_conf

    out_readme = open(circos+".circos.result.legend.txt","w")

    print >>out_readme, "Name for each circle:"
    print >>out_readme, "\t1\tChromosomes"
    circle_cnt = 2
    for i in heatmap_highlights_list:
        print >>out_readme, "\t{}\t{}".format(circle_cnt, i)
        circle_cnt += 1
    #--------------------------------------------------

    if colorD:
        print >>out_readme, "Color legend for highlight:"
        for label, color in colorD.items():
            print >>out_readme, "\t{}\t{}".format(label, color)
    if minv != "ct":
        print >>out_readme, "Color legend for heatmap:"
        print >>out_readme, "Min: {}; Max: {}; \nColor_min: {}; Color_max: {}".format(minv, maxv, 'blue', 'red')
    out_readme.close()

    if os.system(cmd):
        print >>sys.stderr, " Wrong"
        sys.exit(1)
#----------------------------------------

if __name__ == '__main__':
    main()


