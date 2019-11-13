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
    This script is used to generate pymol script for displaying hbonds and listing hbonds.
'''

import sys
import os
from optparse import OptionParser as OP

debug = 0

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -p '/data/prot.pdb,cyp450' -d '/data/chemical1.pdbqt,chem1;/data/chemical2.pdb,chem2' -o output"
    parser = OP(usage=usages)
    parser.add_option("-p", "--protein", dest="prot",
        metavar="FILEIN", help="Supply the absolute path for protein pdb file and protein names used to show in PyMOL. In format like <pdb_path,prot_name>. The shorter name,  the better. No chinese words or spaces allowed.")
    parser.add_option("-d", "--docking-result", dest="docking",
        help="Supply the absolute path for docking result pdbqt file or ligand pdb file and ligand names used to show in PyMOL. In format like <pdbqt_path,ligand_name;pdb_path, ligand_name2;>. The shorter name, the better. No chinese words or spaces allowed.")
    parser.add_option("-o", "--output-prefix", dest="output",
        help="Outpur file prefix.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.prot != None, "A protein/receptor file needed for -i"
    assert options.output != None, "A string needed for -o"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    prot_path, prot_name = options.prot.replace(' ', '').split(',')
    ligandL = [ligand.split(',') for ligand in 
            options.docking.replace(' ', '').split(';') if ligand]
    output_p = options.output
    h_bond_d = open(output_p+'.hbonds_display.pml', 'w')
    h_bond_l = open(output_p+'.hbonds_list.py', 'w') 
    verbose = options.verbose
    global debug
    debug = options.debug
    colorList = [("br9", "br0"), ("brightorange", "carbon"), ("purple", "purpleblue"), ("limon", "limegreen"), ("hotpink", "lightorange")]
    #-----------------------------------
    aDict = {}
    aDict["prot_path"] = prot_path
    aDict['prot_name'] = prot_name
    
    print >>h_bond_d, '''
# The following 4 lines: 
	# 1. load protein structure and rename it
	# 2. add hydrogen (`h_add` uses a primitive algorithm to add hydrogens onto a molecule.)
	# 3. hide protein display
	# 4. show cartoon display for protein
load {d[prot_path]}, {d[prot_name]}
h_add {d[prot_name]}
hide everything, {d[prot_name]}
show cartoon, {d[prot_name]}
cmd.spectrum("count", selection="{d[prot_name]}", byres=1)
'''.format(d=aDict)
    
    print >>h_bond_d, '''#Define a dict mapping protein 3 letters to 1 letter
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', \
'GLY':'G', 'PRO':'P', 'CYS':'C'}    
'''

    for ligand_path, ligand_name in ligandL:
        aDict['ligand_path'] = ligand_path
        aDict['ligand_name'] = ligand_name

        print >>h_bond_d, '''
# The following 6 lines: 
	# 1. load ligand structure and rename it
	# 2. add hydrogen
	# 3. hide ligand display
	# 4. show ligand in sticks mode
	# 5. Set width of stick to 0.15
	# 6. Set atom color: C-white N-blue O-red
load {d[ligand_path]}, {d[ligand_name]}
h_add {d[ligand_name]}
hide everything, {d[ligand_name]}
show sticks, {d[ligand_name]}
set stick_radius, 0.15
util.cbaw {d[ligand_name]}

# The following 2 lines:
	# 1. Set hydrogen donator
	# 2. Set hydrogen accrptor 
	# `select` creates a named selection from an atom selection. 
	# `select name, (selection)`
select h_donator,  (elem n,o and (neighbor hydro))
select h_acceptor, (elem o or (elem n and not (neighbor hydro)))

# The following 4 lines:
	# 1. Create link between ligand_h_acceptor and prot_h_donator  within given distance 3.2
	# 2. Create link between ligand_h_donator  and prot_h_acceptor within given distance 3.2
	#    Set filter 3.6 for ideal geometry and filter 3.2 for minimally acceptable geometry
	# 3. Set red color for ligand_h_acceptor and prot_h_donator 
	# 4. Set blue color for ligand_h_donator  and prot_h_acceptor
	# `distance` creates a new distance object between two selections. It will display all distances within the cutoff. Distance is also used to make hydrogen bonds like `distance hbonds, all, all, 3.2, mode=2`.
	# distance [ name [, selection1 [, selection2 [, cutoff [, mode ]]]]]
distance P2L_{d[ligand_name]}, ({d[ligand_name]} and h_acceptor), ({d[prot_name]} and h_donator), 3.2
distance L2P_{d[ligand_name]}, ({d[ligand_name]} and h_donator), ({d[prot_name]} and h_acceptor), 3.2
color red, P2L_{d[ligand_name]}
color blue, L2P_{d[ligand_name]}
set transparency, 0.4, P2L_{d[ligand_name]}
set transparency, 0.4, L2P_{d[ligand_name]}

# The following 3 lines:
	# 1. Select non-hydro atoms of ligands
	# 2. Select protein atoms within 5A of selected atoms in last step
	# 3. Label alpha-c(ca) of selected residues with residue name and residue position
select sele, {d[ligand_name]} & not hydro

list = []
file = '{d[ligand_path]}_{d[prot_name]}.neighbor_aa.txt'
file_fh = open(file, 'w')
print >>file_fh, "{d[prot_name]} {d[ligand_name]}"
print >>file_fh, "Res_position;Res_name;Res_name;distance"
select seletmp, byres (sele expand 1) & {d[prot_name]}
iterate (name ca & seletmp), list.append((resi, resn, one_letter[resn], '1'))
print >>file_fh, '\\n'.join([';'.join(i) for i in list])
list = []
select seletmp, byres (sele expand 2) & {d[prot_name]}
iterate (name ca & seletmp), list.append((resi, resn, one_letter[resn], '2'))
print >>file_fh, '\\n'.join([';'.join(i) for i in list])
list = []
select seletmp, byres (sele expand 3) & {d[prot_name]}
iterate (name ca & seletmp), list.append((resi, resn, one_letter[resn], '3'))
print >>file_fh, '\\n'.join([';'.join(i) for i in list])
list = []
select seletmp, byres (sele expand 4) & {d[prot_name]}
iterate (name ca & seletmp), list.append((resi, resn, one_letter[resn], '4'))
print >>file_fh, '\\n'.join([';'.join(i) for i in list])
list = []

select sele, byres (sele expand 5) & {d[prot_name]}
label name ca & sele, "%s-%s" % (one_letter[resn],resi)
iterate (name ca & sele), list.append((resi, resn, one_letter[resn], '5')) 
print >>file_fh, '\\n'.join([';'.join(i) for i in list])
file_fh.close()

# The follwing 5 lines
	# 1. Comment out this line
	# 2. Create an object `surrounding_res` to represent selected protein atoms
	#    `create`: creates a new molecule object from a selection. It can also be used to create states in an   existing object. 
	#    `create name, (selection)`
	# 3. Display created surface
    # 4. Set color for surrounding_res: srd_{d[ligand_name]}
	# 5. Set transparency for surrounding_res
	#    Transparency is used to adjust the transparency of Surfaces and Slices.    
	#    `set transparency, F, selection`
#show surface, {d[prot_name]}
create srd_{d[ligand_name]}, sele
show surface, srd_{d[ligand_name]}
show cartoon, srd_{d[ligand_name]}
color grey80, srd_{d[ligand_name]}
set transparency, 0.5, srd_{d[ligand_name]}
'''.format(d=aDict)

    #--------END all ligand-----------
    print >>h_bond_d, '''
# The following 3 lines
	# 1. Set background white
	# 2. Set label color back
	# 3. Hidden hydrogenes
bg white
set label_color, black
hide (hydro)
'''

    h_bond_d.close()

    print >>h_bond_l, '''
from pymol import cmd

def list_hb(selection,selection2=None,cutoff=3.2,angle=55,mode=1,hb_list_name='hbonds'):
    """
    USAGE
  
    list_hb selection, [selection2 (default=None)], [cutoff (default=3.2)],
                       [angle (default=55)], [mode (default=1)],
                       [hb_list_name (default='hbonds')]
  
    The script automatically adds a requirement that atoms in the
    selection (and selection2 if used) must be either of the elements N or
    O.
  
    If mode is set to 0 instead of the default value 1, then no angle
    cutoff is used, otherwise the angle cutoff is used and defaults to 55
    degrees.
  
    e.g.
    To get a list of all H-bonds within chain A of an object
      list_hb 1abc & c. a &! r. hoh, cutoff=3.2, hb_list_name=abc-hbonds
  
    To get a list of H-bonds between chain B and everything else:
      list_hb 1tl9 & c. b, 1tl9 &! c. b
  
    """
    cutoff=float(cutoff)
    angle=float(angle)
    mode=float(mode)
    # ensure only N and O atoms are in the selection
    selection = selection + " & e. n+o"
    if not selection2:
        hb = cmd.find_pairs(selection,selection,mode=mode,cutoff=cutoff,angle=angle)
    else:
        selection2 = selection2 + " & e. n+o"
        hb = cmd.find_pairs(selection,selection2,mode=mode,cutoff=cutoff,angle=angle)
  
    # sort the list for easier reading
    hb.sort(lambda x,y:(cmp(x[0][1],y[0][1])))
  
    for pairs in hb:
        cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'print "\t%1s/%3s`%s/%-4s " % (chain,resn,resi,name),')
        cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'print "%1s/%3s`%s/%-4s " % (chain,resn,resi,name),')
        print "%.2f" % cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))

#cmd.extend("list_hb",list_hb)
#if __name__ == "__main__":
########Default, do not reload proteins and ligands.
########If you do not have these proteins or ligands loaded, 
########Please uncomment lines beginning with 4 well symbols(<####>).
####cmd.load("{d[prot_path]}", "{d[prot_name]}")
cmd.h_add("({d[prot_name]})")
'''.format(d=aDict)

    for ligand_path, ligand_name in ligandL:
        aDict['ligand_path'] = ligand_path
        aDict['ligand_name'] = ligand_name

        print >>h_bond_l, '''
####cmd.load("{d[ligand_path]}","{d[ligand_name]}")
cmd.h_add("({d[ligand_name]})")
h_donator  = "elem n,o & (neighbor hydro)"
h_acceptor = "elem o | (elem n & !(neighbor hydro))"

lacc = "{d[ligand_name]} & (elem o | (elem n & !(neighbor hydro)))"
ldon = "{d[ligand_name]} & (elem n,o & (neighbor hydro))"
pacc = "{d[prot_name]} & (elem o | (elem n & !(neighbor hydro)))"
pdon = "{d[prot_name]} & (elem n,o & (neighbor hydro))"

print "{d[ligand_name]}_2_{d[prot_name]}_hbonds\\n"
list_hb(ldon, pacc, hb_list_name="{d[ligand_name]}_2_{d[prot_name]}_hbonds", mode=0)
print "\\n{d[prot_name]}_2_{d[ligand_name]}_hbonds\\n"
list_hb(lacc, pdon, hb_list_name="{d[prot_name]}_2_{d[ligand_name]}_hbonds", mode=0)
print '\\n-----------------------\\n'
'''.format(d=aDict)

    h_bond_l.close()


if __name__ == '__main__':
    main()
