#!/bin/bash

#set -x
set -e
set -u

usage()
{
cat <<EOF >&2
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to create configuration file for annocript.pl.

${txtbld}OPTIONS${txtrst}:
	-d	Annoscript working directory 
		${bldred}[Relative path, NECESSARY]${txtrst}
	-D	Annoscript intsall directory 
		${bldred}[Default: /MPATHB/soft/Annoscript, NECESSARY]${txtrst}
	-f	Fastaseqs file. 
		Allowed characters [A-za-z0-9\_\-]. Allowed extensions (fa|fasta). 
		Please use a dot only to separate the extension!	
		${bldred}[NECESSARY]${txtrst}
EOF
}

anno_work_dir=
anno_install_dir=/MPATHB/soft/Annoscript
file=

while getopts "hd:D:f:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		d)
			anno_work_dir=$OPTARG
			;;
		D)
			anno_install_dir=$OPTARG
			;;
		f)
			file=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if test -z "${anno_work_dir}"; then
	usage
	exit 1
fi

cat <<END >${anno_work_dir}/config_user.txt
# File configuration for the user of Annocript
# READ CAREFULLY!!!
# This file has been written with a specific sintax. 
# The variables MUST stay in the format: variable = value
# A series of hashes (#########) closes the parameters to read
# Parameters of BLAST programs without a value assigned will not be used (i.e. word_sizeX = )
# When you want to execute something you have to write YES (in upper case) or NO otherwise
# other strings will give error.
##############################

#Allowed characters [A-za-z0-9\_\-]. Allowed extensions (fa|fasta). Please use a dot only to separate the extension!
fastaSeqs = ${file}

#organisms to blast ('all' means all the organisms in UniProt are taken)
#please use 'all' or a file name with organisms names
#Such file must be placed in your working directory (i.e. ann_works)
#Selection of the organisms works only if the TrEMBL database is used!
blastedOrganism = all

#How to extract GO terms: you can choose to extract 
#for proteins ('proteins'), domains ('domains') or for both ('both')
goTermsAss = both

#Steps to perform
doDbCreation = YES
doExecutePrograms = YES
doBuildOutput = YES
extractStatistics = YES

#Analyses to execute
doBlastxSP = YES
doBlastxTRorUf = YES
doRpstblastn = YES
doBlastn = YES
doPortrait = YES
doDna2Pep = YES

#Generation of a GFF database 
useGFFDB = NO
#Write YES, if you want GFF output files. Using NO increases the speed.)
printGFFOutput = NO

#BLASTX and BLASTP PARAMETERS (we use word_size 4 and threshold 18 to reduce computational time)
#(outfmt can be only 0 with this version of Annocript)
#Currently you can only use these parameters. Please ask in the forum if you need others.
word_sizeX = 4
evalueX = 1E-5
num_descriptionsX = 5
num_alignmentsX = 5
max_target_seqsX = 
num_threadsX = 15
thresholdX = 18
matrixX = 

#BLASTN PARAMETERS
word_sizeN = 
evalueN = 0.00001
num_descriptionsN = 1
num_alignmentsN = 1
max_target_seqsN =
num_threadsN = 4
thresholdN =

#RPSBLAST and RPSTBLASTN PARAMETERS 
word_sizeRPS = 
evalueRPS = 0.00001
num_descriptionsRPS = 20
num_alignmentsRPS = 20
max_target_seqsRPS =
thresholdRPS =

#Number of threads for parallel executions (Used only for RPSBLAST)
threads4Parallel = 30

#BLAST results with evalue lower than evalMax will be shown in the tabular output
evalMax = 0.00001

#DNA2PEP PARAMETERS
d2pMode = none

#PLOTS
#Number of top scored elements to show in the plots (maximum is 50)
topToShow = 20
#Type of plot to show [currently you can use only barplot]
plotType = barplot

#Thresholds to be non-coding. They guide the heuristic in Annocript
#Minimum Portrait score
NCThresh = 0.95
#Maximum length of the ORF
NCORFLength = 100
#Minimum length of the transcript
NCSeqLength = 200

#FIXED PARAMETERS (You should set only once)#

#Database account info
mySqlUser = annoscript
mySqlPass = annoscript123

#UNIPROT informations for access
uniprotWebUser = anonymous
uniprotWebPass = tchen@genetics.ac.cn

#Programs Paths
blastPath = /MPATHB/soft/ncbi-blast-2.2.31+/bin/
portraitPath = /MPATHB/soft/Annoscript/DL_PROGRAMS/portrait-1.1/portrait-1.1.pl
dna2pepPath = /MPATHB/soft/Annoscript/DL_PROGRAMS/dna2pep-1.1/dna2pep.py
##############################
END

cat <<END >${anno_work_dir}/folders.txt
`pwd`/${anno_work_dir} ${anno_install_dir}
END
