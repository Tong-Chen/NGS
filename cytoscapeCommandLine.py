#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to generate cytoscape command line and related xml
    file for plotting gene co-expression network.

<network file>

miRNA	Gene	Cor	Pvalue	Padj
rno-miR-1224	Agap2	-0.962840850717	0.0020455488424	1
rno-miR-1224	Cdk5rap2	0.990478240426	0.000135564218109	1
rno-miR-1224	Cntn2	-0.861628416906	0.0273953632623	1
rno-miR-1224	Dag1	0.991534362451	0.000107197175154	1
rno-miR-1224	Dagla	-0.955499849403	0.00292633409495	1
rno-miR-1224	Il15ra	0.23159478558	0.658818747347	1
rno-miR-1224	Map1b	0.981858220501	0.000490700796083	1
rno-miR-1224	Megf8	-0.984493747337	0.000358801615742	1
rno-miR-1224	Nefm	0.986451318654	0.000274106605508	1

<attribute file>

ID	log(ground/space)	Type
rno-miR-1224	-3.03765	miRNA
Agap2	0.658828	mRNA
Cdk5rap2	-1.82057	mRNA
Cntn2	0.290083	mRNA
Dag1	-0.125188	mRNA


<layout>

  layout apply preferred  Execute the preferred layout on a network
  layout attribute-circle  Execute the Attribute Circle Layout on a network
  layout attributes-layout  Execute the Group Attributes Layout on a network
  layout circular  Execute the Circular Layout on a network
  layout cose  Execute the Compound Spring Embedder (CoSE) on a network
  layout degree-circle  Execute the Degree Sorted Circle Layout on a network
  layout force-directed  Execute the Prefuse Force Directed Layout on a network
  layout force-directed-cl  Execute the Prefuse Force Directed OpenCL Layout on a network
  layout fruchterman-rheingold  Execute the Edge-weighted Force directed (BioLayout) on a network
  layout get preferred  Return the current preferred layout
  layout grid  Execute the Grid Layout on a network
  layout hierarchical  Execute the Hierarchical Layout on a network
  layout isom  Execute the Inverted Self-Organizing Map Layout on a network
  layout kamada-kawai  Execute the Edge-weighted Spring Embedded Layout on a network
  layout set preferred  Set the preferred layout
  layout stacked-node-layout  Execute the Stacked Node Layout on a network

<link>

<http://wiki.cytoscape.org/Cytoscape_3/UserManual/Command_Tool>

'''

import sys
import os
from json import dump as json_dump
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool
import numpy as np
import pandas as pd

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

debug = 0


def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="network",
        metavar="FILEIN", help="Network file with first line as header line")
    parser.add_option("-s", "--source-column", dest="source",
        type="int", default=1, help="Default 1 meaning first column of network file.")
    parser.add_option("-t", "--target-column", dest="target",
        type="int", default=2, help="Default 2 meaning second column of network file.")
    parser.add_option("-e", "--edge-attr-column", dest="edge_attr",
        type="int", default=3, help="Default 3 meaning third column of network file. Edge will be colored according to this line.")
    parser.add_option("-I", "--interaction-type-column", dest="inter_type",
        type="int", help="[Optional] No default. Specify a number <n> to indicate the <n>th column of network file.")
    parser.add_option("-l", "--network-layout", dest="layout",
        default="apply preferred", 
        help="Accept <apply preferred (default)>, <force-directed>, <kamada-kawai (spring embed)>, <attribute-circle>, <cose>. See above for more.")
    parser.add_option("-a", "--attribute-file", dest="attribute",
        help="Attribute file with first line as header line")
    parser.add_option("-c", "--node-color-column", dest="node_color",
        type="int", default=2, help="Default 2 meaning the second coulmn of attribute file.")
    parser.add_option("-S", "--node-shape-column", dest="node_shape",
        type="int", default=3, help="Default 3 meaning the third coulmn of attribute file.")
    parser.add_option("-A", "--node-shape-attribute", dest="node_shape_attr",
            help="Given strings in format like <a:DIAMOND;b:ELLIPSE;c:TRIANGLE;d:VEE;e:ROUND_RECTANGLE;f:HEXAGON;g:OCTAGON;miRNA:PARALLELOGRAM> to specify node shape for each value in <node_shape_column>.")
    parser.add_option("-w", "--absolute-path", dest="abs_path",
            help="Given a path to tell cytoscape the position of input file. Normmaly used when running this script in linux and running cytoscape in windows. So a string like <D:/cytoscape/> (windows) or </home/ehbio/cytoscape_data/> (linux) would be suitable.")
    parser.add_option("-o", "--output-prefix", dest="op_prefix",
        help="Default ehbio. Accept a string containing only alphabets.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.network != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def generateNetworkCommand(aDict):
    fh = open(aDict['cyjs_cmd_file'], 'w')
    print >>fh, """
# help command to see detailed usage of each command
network import file dataTypeList="{aDict[network_colType]}" indexColumnTargetInteraction={aDict[target]} indexColumnSourceInteraction={aDict[source]} file="{aDict[network_full]}" firstRowAsColumnNames=true startLoadRow=1;

table import file dataTypeList="{aDict[attribute_colType]}" caseSensitiveNetworCollectionKeys=true caseSensitiveNetworkKeys=true file="{aDict[attribute_full]}" firstRowAsColumnNames=true startLoadRow=1 keyColumnIndex=1;

vizmap load file file="{aDict[style_xml_full]}";
vizmap apply styles={aDict[style_name]};

layout {aDict[layout]};

view fit content  

network export options="Cytoscape.js JSON (*.cyjs)" OutputFile="{aDict[output_prefix_full]}"
view export OutputFile="{aDict[output_prefix_full]}" options=PDF

""".format(aDict=aDict)
    fh.close()
#----------------------------------------------

def generateXml(aDict, node_shape_attrD):

    fh = open(aDict['style_xml'], 'w')

    print >>fh, """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<vizmap documentVersion="3.0" id="VizMap-2017_08_20-11_07">
    <visualStyle name="{aDict[style_name]}">
        <network>
            <visualProperty name="NETWORK_DEPTH" default="0.0"/>
            <visualProperty name="NETWORK_EDGE_SELECTION" default="true"/>
            <visualProperty name="NETWORK_HEIGHT" default="400.0"/>
            <visualProperty name="NETWORK_WIDTH" default="550.0"/>
            <visualProperty name="NETWORK_SCALE_FACTOR" default="1.0"/>
            <visualProperty name="NETWORK_TITLE" default=""/>
            <visualProperty name="NETWORK_CENTER_Y_LOCATION" default="0.0"/>
            <visualProperty name="NETWORK_BACKGROUND_PAINT" default="#FFFFFF"/>
            <visualProperty name="NETWORK_SIZE" default="550.0"/>
            <visualProperty name="NETWORK_CENTER_X_LOCATION" default="0.0"/>
            <visualProperty name="NETWORK_NODE_SELECTION" default="true"/>
            <visualProperty name="NETWORK_CENTER_Z_LOCATION" default="0.0"/>
        </network>
        <node>
            <dependency name="nodeCustomGraphicsSizeSync" value="true"/>
            <dependency name="nodeSizeLocked" value="false"/>
            <visualProperty name="NODE_Z_LOCATION" default="0.0"/>
            <visualProperty name="NODE_CUSTOMPAINT_6" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_6, name=Node Custom Paint 6)"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_5" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_9" default="50.0"/>
            <visualProperty name="NODE_CUSTOMPAINT_5" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_5, name=Node Custom Paint 5)"/>
            <visualProperty name="NODE_TRANSPARENCY" default="255"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_8" default="50.0"/>
            <visualProperty name="NODE_CUSTOMPAINT_7" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_7, name=Node Custom Paint 7)"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_3" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_7" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
            <visualProperty name="NODE_LABEL_COLOR" default="#000000"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_9" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_CUSTOMPAINT_9" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_9, name=Node Custom Paint 9)"/>
            <visualProperty name="NODE_LABEL_TRANSPARENCY" default="255"/>
            <visualProperty name="NODE_X_LOCATION" default="0.0"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_8" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
            <visualProperty name="NODE_TOOLTIP" default=""/>
            <visualProperty name="NODE_CUSTOMPAINT_4" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_4, name=Node Custom Paint 4)"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_4" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_7" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_1" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_CUSTOMPAINT_3" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_3, name=Node Custom Paint 3)"/>
            <visualProperty name="NODE_WIDTH" default="90.0"/>
            <visualProperty name="NODE_Y_LOCATION" default="0.0"/>
            <visualProperty name="NODE_LABEL_POSITION" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_CUSTOMPAINT_8" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_8, name=Node Custom Paint 8)"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_2" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_6" default="C,C,c,0.00,0.00"/>
            <visualProperty name="COMPOUND_NODE_PADDING" default="10.0"/>
            <visualProperty name="COMPOUND_NODE_SHAPE" default="ROUND_RECTANGLE"/>
            <visualProperty name="NODE_VISIBLE" default="true"/>
            <visualProperty name="NODE_LABEL_FONT_FACE" default="SansSerif.plain,plain,12"/>
            <visualProperty name="NODE_HEIGHT" default="35.0"/>
            <visualProperty name="NODE_LABEL_FONT_SIZE" default="25"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_3" default="50.0"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_5" default="50.0"/>
            <visualProperty name="NODE_BORDER_STROKE" default="SOLID"/>
            <visualProperty name="NODE_BORDER_PAINT" default="#CCCCCC"/>
            <visualProperty name="NODE_SELECTED_PAINT" default="#FFFF00"/>
            <visualProperty name="NODE_CUSTOMPAINT_2" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_2, name=Node Custom Paint 2)"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_8" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_6" default="50.0"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_3" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_SELECTED" default="false"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_1" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
            <visualProperty name="NODE_CUSTOMPAINT_1" default="DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_1, name=Node Custom Paint 1)"/>
            <visualProperty name="NODE_LABEL_WIDTH" default="200.0"/>
            <visualProperty name="NODE_BORDER_WIDTH" default="0.0"/>
            <visualProperty name="NODE_SHAPE" default="ROUND_RECTANGLE">
                <discreteMapping attributeType="string" attributeName="Type">""".format(aDict=aDict)

    for item, shape in node_shape_attrD.items():
        print >>fh, " "*20+"<discreteMappingEntry value=\""+shape+'" attributeValue="'+item+'"/>'

    print >>fh, """                </discreteMapping>
            </visualProperty>
            <visualProperty name="NODE_NESTED_NETWORK_IMAGE_VISIBLE" default="true"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_4" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
            <visualProperty name="NODE_LABEL" default="">
                <passthroughMapping attributeType="string" attributeName="name"/>
            </visualProperty>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_2" default="50.0"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_6" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_9" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
            <visualProperty name="NODE_FILL_COLOR" default="#89D0F5">
                <continuousMapping attributeType="float" attributeName="{aDict[node_color_name]}">
                    <continuousMappingPoint lesserValue="#000000" greaterValue="#0000FF" equalValue="#0000FF" attrValue="{aDict[node_color_min]}"/>
                    <continuousMappingPoint lesserValue="#FFFFFF" greaterValue="#FFFFFF" equalValue="#FFFFFF" attrValue="0"/>
                    <continuousMappingPoint lesserValue="#FFFF66" greaterValue="#FFFFFF" equalValue="#FFFF66" attrValue="{aDict[node_color_max]}"/>
                </continuousMapping>
            </visualProperty>
            <visualProperty name="NODE_SIZE" default="35.0"/>
            <visualProperty name="NODE_BORDER_TRANSPARENCY" default="255"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_1" default="50.0"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_4" default="50.0"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_POSITION_2" default="C,C,c,0.00,0.00"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_SIZE_7" default="50.0"/>
            <visualProperty name="NODE_PAINT" default="#1E90FF"/>
            <visualProperty name="NODE_DEPTH" default="0.0"/>
            <visualProperty name="NODE_CUSTOMGRAPHICS_5" default="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],"/>
        </node>
        <edge>
            <dependency name="arrowColorMatchesEdge" value="false"/>
            <visualProperty name="EDGE_LINE_TYPE" default="SOLID"/>
            <visualProperty name="EDGE_WIDTH" default="2.0"/>
            <visualProperty name="EDGE_STROKE_UNSELECTED_PAINT" default="#848484">
                <continuousMapping attributeType="float" attributeName="{aDict[edge_attr_name]}">
                    <continuousMappingPoint lesserValue="#000000" greaterValue="#66FF66" equalValue="#66FF66" attrValue="{aDict[edge_attr_min]}"/>
                    <continuousMappingPoint lesserValue="#CCCCFF" greaterValue="#CCCCFF" equalValue="#CCCCFF" attrValue="0"/>
                    <continuousMappingPoint lesserValue="#FF3333" greaterValue="#FFFFFF" equalValue="#FF3333" attrValue="{aDict[edge_attr_max]}"/>
                </continuousMapping>
            </visualProperty>
            <visualProperty name="EDGE_LABEL" default=""/>
            <visualProperty name="EDGE_LABEL_TRANSPARENCY" default="255"/>
            <visualProperty name="EDGE_VISIBLE" default="true"/>
            <visualProperty name="EDGE_LABEL_FONT_SIZE" default="10"/>
            <visualProperty name="EDGE_LABEL_COLOR" default="#000000"/>
            <visualProperty name="EDGE_SELECTED_PAINT" default="#FF0000"/>
            <visualProperty name="EDGE_CURVED" default="true"/>
            <visualProperty name="EDGE_STROKE_SELECTED_PAINT" default="#FF0000"/>
            <visualProperty name="EDGE_PAINT" default="#323232"/>
            <visualProperty name="EDGE_TRANSPARENCY" default="255"/>
            <visualProperty name="EDGE_TARGET_ARROW_UNSELECTED_PAINT" default="#000000"/>
            <visualProperty name="EDGE_BEND" default=""/>
            <visualProperty name="EDGE_TARGET_ARROW_SHAPE" default="NONE"/>
            <visualProperty name="EDGE_UNSELECTED_PAINT" default="#404040"/>
            <visualProperty name="EDGE_SOURCE_ARROW_UNSELECTED_PAINT" default="#000000"/>
            <visualProperty name="EDGE_LABEL_WIDTH" default="200.0"/>
            <visualProperty name="EDGE_LABEL_FONT_FACE" default="Dialog.plain,plain,10"/>
            <visualProperty name="EDGE_SOURCE_ARROW_SHAPE" default="NONE"/>
            <visualProperty name="EDGE_SELECTED" default="false"/>
            <visualProperty name="EDGE_TARGET_ARROW_SELECTED_PAINT" default="#FFFF00"/>
            <visualProperty name="EDGE_SOURCE_ARROW_SELECTED_PAINT" default="#FFFF00"/>
            <visualProperty name="EDGE_TOOLTIP" default=""/>
        </edge>
    </visualStyle>
</vizmap>
""".format(aDict=aDict)
    fh.close()
#---------------------------------

def readNetworkAttribute(network):
    network_data = pd.read_table(network, sep="\t", header=0)
    dtypes = network_data.dtypes
    col_nameL = list(dtypes.index)

    typeL = []
    for col_name in col_nameL:
        cur_type = dtypes[col_name].type
        if str(cur_type) == "<type 'numpy.datetime64'>":
            col_type = 'string'
        elif issubclass(cur_type, np.datetime64):
            col_type = 'string'
        elif issubclass(cur_type, (np.integer, np.bool_)):
            col_type = 'int'
        elif issubclass(cur_type, np.floating):
            col_type = 'double'
        else:
            col_type = 'string'
        typeL.append(col_type)
    typestr = ','.join(typeL)
    network_data_summary = network_data.describe()
    return col_nameL, typestr, network_data_summary
#-------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    network      = options.network
    source       = options.source
    target       = options.target
    if options.inter_type:
        inter_type   = options.inter_type - 1
    else:
        inter_type = ''
    layout       = options.layout
    edge_attr    = options.edge_attr - 1
    attribute    = options.attribute
    node_color   = options.node_color - 1
    node_shape   = options.node_shape - 1
    node_shape_attr = options.node_shape_attr
    if node_shape_attr:
        node_shape_attrD = dict([i.split(':') for i in node_shape_attr.split(';')])
    else:
        node_shape_attrD = {}
    op_prefix    = options.op_prefix
    abs_path     = options.abs_path
    verbose      = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    aDict = {'source': source, 'target': target, 'inter_type':inter_type, 
            'network': network, 'layout': layout, 
            'network_full': abs_path+'/'+os.path.split(network)[1], 
            'edge_attr': edge_attr, 'attribute': attribute, 
            'attribute_full': abs_path+'/'+os.path.split(attribute)[1], 
            'node_color': node_color, 'node_shape': node_shape, 
            'style_xml': op_prefix+'.style.xml', 
            'style_xml_full': abs_path+'/'+op_prefix+'.style.xml', 
            'style_name': op_prefix, 
            'output_prefix': op_prefix, 
            'output_prefix_full': abs_path+'/'+op_prefix, 
            'cyjs_cmd_file': op_prefix+'.net.cmd'
            }

    #-----------------------------------
    network_colName, network_colType, network_data_summary \
            = readNetworkAttribute(network)
    attribute_colName, attribute_colType, attribute_data_summary \
            = readNetworkAttribute(attribute)
    aDict['network_colType'] = network_colType
    aDict['attribute_colType'] = attribute_colType

    aDict['edge_attr_name'] = network_colName[edge_attr]
    aDict['edge_attr_min'] = network_data_summary[network_colName[edge_attr]][3]
    aDict['edge_attr_max'] = network_data_summary[network_colName[edge_attr]][7]

    aDict['node_color_name'] = attribute_colName[node_color]
    #print >>sys.stderr, aDict
    aDict['node_color_min'] = attribute_data_summary[attribute_colName[node_color]][3]
    aDict['node_color_max'] = attribute_data_summary[attribute_colName[node_color]][7]

    if length(attribute_colName) > 3:
        aDict['node_shape_name'] = attribute_colName[node_shape]
    else:
        aDict['node_shape_name'] = 0


    generateNetworkCommand(aDict)
    generateXml(aDict, node_shape_attrD)
    


    if verbose:
        print >>sys.stderr,\
            "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


