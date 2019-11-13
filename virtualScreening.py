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
    This is designed to do virtual screening between proteins and ligands.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
from numpy import mean, median, std, ptp
from multiprocessing.dummy import Pool as ThreadPool

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
    usages = "%prog -p 'prot1.pdb,prot2.pdbqt' -l 'm1.mol2,m2.pdb,m3.pdbqt'"
    parser = OP(usage=usages)
    parser.add_option("-p", "--prot", dest="prot",
        metavar="FILEIN", help="<,> separated one or multiple \
protein PDB or PDBQT file defined \
by suffix. If PDB files are given, it will be translated into AutoVina \
required PDBQT files.")
    parser.add_option("-l", "--ligand", dest="ligand",
        help="<,> separated one or multiple ligand file in mol2, pdb or pdbqt file. Files not in pdbqt format will be transferred to pdbqt file.")
    parser.add_option("-c", "--center", dest="center",
        default='auto', 
        help="The coordinate of grid center. Default the program will \
take the center of protein as grid center. Both mean value and median \
value will be used as center. One can also supply \
center coordinates (x,y,z) in format like <10,5,6;20,-10,30> \
for each protein. Default auto.")
    parser.add_option("-s", "--size", dest="size",
        default='auto', 
        help="The size of the grid box (Angstroms). \
Default the program will compute the grid size to include \
as many residues. It has three dimentions, center+std, center+2*std, \
center+ptp/2. One can also supply \
size for (x,y,z) in format like <30,30,30;25,30,40> \
for each protein. Default auto.")
    parser.add_option("-n", "--number-of-runs-for-each-pair", 
        default=3, type=int, 
        dest="num_run", help="Specify the number of independent runs to perform to get the best result. Default 3.")
    parser.add_option("-f", "--fine-tune", dest="fineTune",
        default=False, action="store_true", 
        help="Fine tune the docking model (shrink pocket size after one run of docking)")
    parser.add_option("-o", "--output-dir", dest="out_dir",
        help="The path for output directory")
    parser.add_option("-t", "--thread", dest="thread",
        default=1, type="int", help="Specify number of threads running to accelerate program running. Default 1.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.prot != None, "A losted of proteins needed for -p"
    return (options, args)
#--------------------------------------------------------------------

def protPrepare(protL):
    len_prot = len(protL)
    for i in range(len_prot):
        prot = protL[i]
        file_n, suffix = prot.rsplit('.', 1)
        if suffix.lower() != 'pdbqt':
            pdbqt = file_n + '.pdbqt'
            protL[i] = pdbqt
            if not os.path.exists(pdbqt):
                cmd = ['prepare_receptor4.py -r', prot, '-o',  
                        pdbqt, '-A hydrogens']
                os.system(' '.join(cmd))
#-----------------------------------
def ligandPrepare(ligandL):
    len_ligand = len(ligandL)
    for i in range(len_ligand):
        ligand = ligandL[i]
        file_n, suffix = ligand.rsplit('.', 1)
        if suffix.lower() != 'pdbqt':
            pdbqt = file_n + '.pdbqt'
            ligandL[i] = pdbqt
            if not os.path.exists(pdbqt):
                cmd = ['prepare_ligand4.py -l', ligand, '-o',  
                        pdbqt, '-A hydrogens']
                os.system(' '.join(cmd))
#---------------------------------------------------------
def determineCenterSize(protL):
    '''
    aDict = {'prot': 
              [[(mean_x, mean_y, mean_z), (median_x, median_y, median_z)], 
               [(std_x,  std_y,  std_z),  (ptp_x,    ptp_y,    ptp_z   )]
              ]
            }
    '''
    aDict = {}
    for prot in protL:
        x = []
        y = []
        z = []
        for line in open(prot):
            if line.startswith('ATOM'):
                x.append(float(line[30:38]))
                y.append(float(line[38:46]))
                z.append(float(line[46:54]))
        #-------------------------------------------------
        std_x = std(x)
        std_y = std(y)
        std_z = std(z)
        aDict[prot] = [[(str(mean(x)), str(mean(y)), str(mean(z))), 
                        (str(median(x)), str(median(y)), str(median(z)))], 
                       [(str(std_x), str(std_y), str(std_z)), 
                        (str(std_x*2), str(std_y*2), str(std_z*2)), 
                        (str(ptp(x)/2), str(ptp(y)/2), str(ptp(z)/2))]
                      ]
        if debug:
            print >>sys.stderr, aDict[prot]
    return aDict
#---------------------------------------------------

def readLog(count, logfile):
    '''
    logL = [[count, mode, energy], [count, mode, energy] ]
    '''
    logL = []
    fh = open(logfile)
    line = fh.readline()
    while not line.startswith('---'):
        line = fh.readline()

    line = fh.readline()

    while not line.startswith("Writing output"):
        lineL = line.strip().split()
        #print >>sys.stderr, lineL
        assert len(lineL) == 4, "Wrong log for {}".format(logfile)
        mode = int(lineL[0])
        energy = float(lineL[1])
        logL.append([count, mode, energy])
        line = fh.readline()
    return logL
#----------------------------------
def readDockingResult(countD, resultF, coord=0):
    '''
    countD = {mode1: pdbqt, mode2:pdbqt}
    pdbqtCoordD = {mode1:[x, y, z], mode2:[x, y, z]} 
    '''
    pdbqtCoordD = {}
    for line in open(resultF):
        line = line.rstrip()
        if line.startswith('MODEL '):
            mode = int(line.split()[1])
            tmpL = []
            countD[mode] = tmpL
            if coord:
                x = []
                y = []
                z = []
                pdbqtCoordD[mode] = [x,y,z]
        else:
            tmpL.append(line)
            if coord and (line.startswith('HETATM') 
                    or line.startswith('ATOM')):
                x.append(float(line[30:38]))
                y.append(float(line[38:46]))
                z.append(float(line[46:54]))
    return pdbqtCoordD
                
#----------------------------------
def extractDockingResult(file, modeL):
    pass

def fineTune(prot_name, ligand_name, runD):
    '''
    runD = {
        1: [(center_x,y,z),(size_x,y,z),log,result], 
        2: [(center_x,y,z),(size_x,y,z),log,result]
    }

    logL = [[run, mode, energy], [run, mode, energy] ]

    resultD = {run:{mode1: [x, y, z], mode2:[x, y, z]}}
    '''
    newX = []
    newY = []
    newZ = []
    logL = []
    resultD = {}
    for run, valueL in runD.items():
        logL.extend(readLog(run, valueL[2]))
        tmpD = {}
        pdbqtCoordD = readDockingResult(tmpD, valueL[3], 1)
        resultD[run] = pdbqtCoordD
        tmpD = {}
    #----------------------------------
    logL.sort(key=lambda x: x[2])
    
    cpd_ptp_x = cpd_ptp_y = cpd_ptp_z = 0

    max = logL[0][2]
    for log in logL:
        run, mode, energy = log
        if energy - max > 2:
            break
        x, y, z = resultD[run][mode]   
        ptp_x,ptp_y,ptp_z = ptp(x),ptp(y),ptp(z)
        cpd_ptp_x = ptp_x if cpd_ptp_x < ptp_x else cpd_ptp_x
        cpd_ptp_y = ptp_y if cpd_ptp_y < ptp_y else cpd_ptp_y
        cpd_ptp_z = ptp_z if cpd_ptp_z < ptp_z else cpd_ptp_z

        if debug:
            print >>sys.stderr, "RUN {} MODE {}\n\t{}\n\t{}\n\t{}".format(run,mode,x,y,z)
        newX.extend(x)
        newY.extend(y)
        newZ.extend(z)
    #--------------------------------------
    x, y, z = newX, newY, newZ
    if debug:
        print >>sys.stderr, "x {}\ny {}\nz {}".format(x, y, z)
    ptp_x,ptp_y,ptp_z = ptp(x),ptp(y),ptp(z)
    std_x, std_y, std_z = std(x), std(y), std(z)
    return [[(str(mean(x)), str(mean(y)), str(mean(z))), 
             (str(median(x)), str(median(y)), str(median(z)))], 
            [(str(cpd_ptp_x), str(cpd_ptp_y), str(cpd_ptp_z)), 
             (str(cpd_ptp_x+std_x), str(cpd_ptp_y+std_y), str(cpd_ptp_z+std_z)), 
             (str(ptp_x), str(ptp_y), str(ptp_z)), 
             (str(ptp_x+std_x), str(ptp_y+std_y), str(ptp_z+std_z)), 
            ]
           ]
#--------------------------------------------

def parseResult(prot_name, ligand_name, out_dir, runD):
    '''
    runD = {
        1: [(center_x,y,z),(size_x,y,z),log,result], 
        2: [(center_x,y,z),(size_x,y,z),log,result]
    }

    logL = [[run, mode, energy], [run, mode, energy] ]

    resultD = {run:{mode1: pdbqt, mode2:pdbqt}}
    '''
    logL = []
    resultD = {}
    for run, valueL in runD.items():
        logL.extend(readLog(run, valueL[2]))
        tmpD = {}
        readDockingResult(tmpD, valueL[3])
        resultD[run] = tmpD
    #----------------------------------
    logL.sort(key=lambda x: x[2])
    logFile = out_dir + '/' + prot_name + '_' + ligand_name + '.log'
    pdbqtFile = out_dir + '/' + prot_name + '_' + ligand_name + '.pdbqt'
    
    logFH = open(logFile, 'w')
    pdbqtFH = open(pdbqtFile, 'w')

    print >>logFH, "MODE\taffinity ENERGY(kcal/mol)\tcenter\tsize\tOriginal_mode"
    mode = 1
    max = logL[0][2]
    for log in logL:
        run = log[0]
        energy = log[2]
        if energy - max > 3:
            break
        #run_p: running parameter
        run_p = runD[run]
        center = run_p[0]
        size   = run_p[1]
        print >>logFH, '{}\t{}\t{}\t{}\t{}_{}'.format(mode,energy,center,size, run, log[1])
        print >>pdbqtFH, 'MODEL {}'.format(mode)
        print >>pdbqtFH, '\n'.join(resultD[run][log[1]])
        mode += 1
    logFH.close()
    pdbqtFH.close()
#--------------------------------------------
def vina_one(prot, prot_name, ligand, datadir, centerL, sizeL, runD, count, num_run):
    cmdL = []
    for center in centerL:
        for size in sizeL:
            for repeat in range(num_run):
                log = datadir+str(count)+'.log'
                dock_output = datadir+str(count)+'.pdbqt'
                runD[count] = [center, size, log, dock_output]
                cmd = ["vina --receptor", prot, 
                        "--ligand", ligand, 
                        "--center_x", center[0], 
                        "--center_y", center[1], 
                        "--center_z", center[2], 
                        "--size_x", size[0], 
                        "--size_y", size[1], 
                        "--size_z", size[2], 
                        "--out", dock_output, 
                        "--log", log, 
                        "--cpu 8",
                        "--seed", str(count), 
                        "--exhaustiveness 20", 
                        "--num_modes 50",
                        ]
                cmd = ' '.join(cmd)
                if debug:
                    cmdL.append(cmd)
                if not os.path.exists(dock_output):
                    #print >>sys.stderr, ">>Running {}".format(cmd)
                    #os.system(cmd)
                    cmdL.append(cmd)
                else:
                    print >>sys.stderr, "{} already exists".format(dock_output)
                count += 1
            #----END repeat-------------
        #--------END size---------------------------
    #-----------------------------------------------------------
    return count,cmdL
#-----------END vina_one------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    protL = options.prot.split(',')
    ligandL = options.ligand.split(',')
    out_dir = options.out_dir
    num_run = options.num_run
    thread  = options.thread
    center = options.center
    DofineTune = options.fineTune
    if center == "auto":
        centerL = []
    else:
        centerL = [i.split(',') for i in center.split(';')]

    size = options.size
    if size == "auto":
        sizeL = []
    else:
        sizeL = [i.split(',') for i in size.split(';')]

    verbose = options.verbose
    global debug
    debug = options.debug

    protPrepare(protL)
    ligandPrepare(ligandL)


    '''
    coordDict = {'prot': 
      #center [[(mean_x, mean_y, mean_z), (median_x, median_y, median_z)], 
      #size    [(std_x,  std_y,  std_z),  (ptp_x,    ptp_y,    ptp_z   )]
              ]
            }
    '''
    coordDict = determineCenterSize(protL)
    
    '''
    recordDict = {prot: {
                   ligand :{
                    'run1': [(center_x,y,z),(size_x,y,z),log,result]
                   }}
                }
    '''
    recordDict = {}

    for prot in protL:
        prot_name = os.path.split(prot)[1][:-6]
        if not centerL:
            centerL, sizeL = coordDict[prot]
        print >>sys.stderr, "Center coordinates: {} and size: {}".format(centerL, sizeL)
        ligandD = {}
        recordDict[prot_name] = ligandD
        for ligand in ligandL:
            ligand_name = os.path.split(ligand)[1][:-6]
            datadir = out_dir + '/' + prot_name+'__'+ligand_name+'/'
            os.system("mkdir -p "+datadir)
            runD = {}
            ligandD[ligand_name] = runD
            count = 1
            count,cmdL = vina_one(prot, prot_name, ligand, datadir, 
                    centerL, sizeL, runD, count, num_run)
            if debug:
                print >>sys.stderr, centerL
                print >>sys.stderr, sizeL
                print >>sys.stderr, '\n'.join(cmdL)
            else:
                if thread > 1:
                    pool = ThreadPool(thread) # thread represents thread_num
                    result = pool.map(os.system, cmdL)
                    pool.close()
                    pool.join()
                else:
                    for cmd in cmdL:
                        os.system(cmd)
            if DofineTune:
                centerL, sizeL = fineTune(prot_name, ligand_name, runD)
                print >>sys.stderr, "New center coordinates: {} and size: {}".format(centerL, sizeL)
                count,cmdL = vina_one(prot, prot_name, ligand, datadir, 
                        centerL, sizeL, runD, count, num_run)
                if debug:
                    print >>sys.stderr, centerL
                    print >>sys.stderr, sizeL
                    print >>sys.stderr, '\n'.join(cmdL)
                else:
                    if thread > 1:
                        pool = ThreadPool(thread) # thread represents thread_num
                        result = pool.map(os.system, cmdL)
                        pool.close()
                        pool.join()
                    else:
                        for cmd in cmdL:
                            os.system(cmd)
                    #----------------------------------
                #----------------------------------
            #---------------------------------------------
#            for center in centerL:
#                for size in sizeL:
#                    for repeat in range(num_run):
#                        log = datadir+str(count)+'.log'
#                        dock_output = datadir+str(count)+'.pdbqt'
#                        runD[count] = [center, size, log, dock_output]
#                        cmd = ["vina --receptor", prot, 
#                                "--ligand", ligand, 
#                                "--center_x", center[0], 
#                                "--center_y", center[1], 
#                                "--center_z", center[2], 
#                                "--size_x", size[0], 
#                                "--size_y", size[1], 
#                                "--size_z", size[2], 
#                                "--out", dock_output, 
#                                "--log", log, 
#                                "--cpu 10",
#                                "--seed", str(count), 
#                                "--exhaustiveness 20", 
#                                "--num_modes 50",
#                                ]
#                        cmd = ' '.join(cmd)
#                        print >>sys.stderr, cmd
#                        os.system(cmd)
#                        count += 1
#                    #----END repeat-------------
#                #--------END size---------------------------
#            #-----------------------------------------------------------
            parseResult(prot_name, ligand_name, out_dir, runD)
    ###--------multi-process------------------
    ###--------multi-process------------------
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


