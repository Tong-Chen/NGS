#!/home/ct/anaconda3/envs/py3.5/bin/python3
# -*- coding: utf-8 -*-
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:

'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

from wordcloud import WordCloud, STOPWORDS
import codecs
import jieba
#import jieba.analyse as analyse
from scipy.misc import imread
from os import path
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print(desc)
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="INPUT text file")
    parser.add_option("-m", "--mask-pic", dest="mask_pic",
        help="A picture used for shaping word cloud. Optional.")
    parser.add_option("-s", "--stopwords", dest="stopwords",
        default="/MPATHB/self/NGS/stopwords", help="Default stopwords in ")
    parser.add_option("-f", "--font", dest="font",
        default="/MPATHB/self/NGS/FZ_School.ttf", help="Fonts (default /MPATHB/self/NGS/FZ_School.ttf)")
    parser.add_option("-o", "--output", dest="output",
        help="Output file name with extensions. If not specified, it would in <input_file.png>")
    parser.add_option("-M", "--maxwords", dest="maxwords",
        type = "int", default=500, 
        help="Maximum allowed words for showing. Default 500.")
    parser.add_option("-F", "--max-font-size", dest="max_font_size",
        type = "int", default=100, 
        help="Maximum font size. Default 100.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

# 绘制词云
def draw_wordcloud(txt, output, font, max_words=500, max_font_size=100, mask_pic=None, stopwords=STOPWORDS):
    #读入一个txt文件
    comment_text = open(txt, 'r').read()
    #结巴分词，生成字符串，如果不通过分词，无法直接生成正确的中文词云
    cut_text = " ".join(jieba.cut(comment_text))
    if mask_pic:
        color_mask = imread(mask_pic) # 读取背景图片
    cloud = WordCloud(
        #设置字体，不指定就会出现乱码
        font_path=font,
        #设置背景色
        background_color='white',
        #词云形状
        mask=color_mask,
        #允许最大词汇
        max_words=max_words,
        #最大号字体
        max_font_size=max_font_size, 
        #
        stopwords=stopwords
    )
    word_cloud = cloud.generate(cut_text) # 产生词云
    word_cloud.to_file(output) #保存图片
    #  显示词云图片
    #plt.imshow(word_cloud)
    #plt.axis('off')
    #plt.show()
#---------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    mask_pic = options.mask_pic
    if not mask_pic:
        mask_pic = None
    stopwords = options.stopwords
    maxwords  = options.maxwords
    font = options.font
    output = options.output
    max_font_size = options.max_font_size
    if not output:
        output = file + '.png'
    verbose = options.verbose
    global debug
    debug = options.debug
    global STOPWORDS
    if stopwords:
        STOPWORDS = STOPWORDS.union(set([line.strip() for line in open(stopwords)]))
    print(STOPWORDS)
    #-----------------------------------
    draw_wordcloud(txt=file, output=output, font=font, max_font_size=max_font_size, 
        max_words=maxwords, mask_pic=mask_pic, stopwords=STOPWORDS)
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    ###---------profile the program---------


