#!/bin/bash

set -x
set -e
set -u

usage()
{
cat <<EOF >&2
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to do PASA transcript annotation.

${txtbld}OPTIONS${txtrst}:
	-p	PASA home${bldred}[Default: /MPATHB/soft/PASApipeline-2.0.2/; 
		NECESSARY]${txtrst}
	-d	Database name${bldred}[NECESSARY]${txtrst}
	-m	Max intron length${bldred}[Default 100000, NECESSARY]${txtrst}
	-c	Cufflinks gtf${bldred}Optional${txtrst}
	-g	The FASTA file for genome. Usually in /MPATHB/resource/ucsc.
	-t	Transcripts de novo assembled using Trinity.
	-T	Transcripts genome-guide assemble using Trinity.
	-G	GFF3 formated protein-coding gene annotation file.
		${bldred}1. Only protein-coding genes should be supplied
		2. ${pasa_home}/misc_utilities/pasa_gff3_validator.pl
		input.gff3 to validate the compatibility of gff file.	
		${txtrst}
	-s	Strand-specific data, accept "<unstrand>" and <FR> or <RF>.
	-r	Rerun Launch_PASA_pipeline due to errors. Accept a number to
		specify the re-run index.
EOF
}

file=
pasa_home="/MPATHB/soft/PASApipeline-2.0.2/"
db=
max_intron_len=100000
cuff=''
genome=''
trinity_denovo=''
trinity_genome=''
gff3=''
strand=
rerun_index=

while getopts "hc:d:g:G:m:p:r:s:t:T:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		c)
			cuff=$OPTARG
			;;
		d)
			db=$OPTARG
			;;
		g)
			genome=$OPTARG
			;;
		G)
			gff3=$OPTARG
			;;
		m)
			max_intron_len=$OPTARG
			;;
		p)
			pasa_home=$OPTARG
			;;
		r)
			rerun_index=$OPTARG
			;;
		s)
			strand=$OPTARG
			;;
		t)
			trinity_denovo=$OPTARG
			;;
		T)
			trinity_genome=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z ${trinity_denovo} ]; then
	usage
	exit 1
fi

if [ -n $cuff ]; then
	cuff="--cufflinks_gtf $cuff"
fi

if [ "${strand}" == "unstrand" ]; then
	orient=
elif [ "${strand}" == "RF" ] || [ ${strand} == "FR" ];then
	orient="--transcribed_is_aligned_orient"
else
	echo "Unknown parameter: -S ${strand}" >&2
	exit 1
fi

if [ "${rerun_index}" != "" ]; then
	echo "Rerun PASA alignment assembl --`date` " >&2
	${pasa_home}/scripts/Launch_PASA_pipeline.pl \
		-c pasa/pasa.alignAssembly.Template.txt \
		-R --ALIGNERS "gmap,blat" -N 1 \
		--MAX_INTRON_LENGTH ${max_intron_len} ${cuff} --CPU 10 \
		-g pasa/genome.fa -t pasa/transcripts.fasta --TDN pasa/trinity_denovo.accs \
		${orient} -L --annots_gff3 ${gff3} \
		--gene_overlap 50 -s ${rerun_index} -d
else
	echo "Begin $0 $@ -- `date`" >&2
	echo " " >&2

	# Concatenate the Trinity.fasta and Trinity.GG.fasta files
	# into a single transcripts.fasta file
	cat ${trinity_denovo} ${trinity_genome} >pasa/transcripts.fasta

	# Create a file containing the list of transcript accessions that
	# correspond to the Trinity de novo assembly (full de novo,  not
	# genome-guided)
	${pasa_home}/misc_utilities/accession_extractor.pl <${trinity_denovo} >pasa/trinity_denovo.accs

	# SET alignment template
	/bin/cp -f ${pasa_home}/pasa_conf/pasa.alignAssembly.Template.txt pasa
	sed -i -e "s/<__MYSQLDB__>/${db}/" -e 's/<__MIN_PERCENT_ALIGNED__>/30/' \
		-e 's/<__MIN_AVG_PER_ID__>/95/' pasa/pasa.alignAssembly.Template.txt

	# Generate mysql username and passwd
	cat <<END >pasa/init.sql
	drop database if exists ${db};
	#create database ${db};
	grant all privileges on ${db}.* to 'pasa'@'localhost' identified by 'pasa_annotate123';
	flush privileges;

END

	mysql -uroot <pasa/init.sql

	##---------Unused--------------------------------
	## seqclean
	#seqclean pasa/transcripts.fasta
	##---------Unused--------------------------------

	##---------Unused--------------------------------
	## Load pre-existing protein-coding gene annotation
	#
	#${pasa_home}/scripts/Load_Current_Gene_Annotations.dbi \
	#	-c pasa/pasa.alignAssembly.Template -g ${genome} -p ${gff3}
	##---------Unused--------------------------------
	ln -sf ${genome} pasa/genome.fa

	echo "PASA alignment assembl --`date` " >&2
	${pasa_home}/scripts/Launch_PASA_pipeline.pl \
		-c pasa/pasa.alignAssembly.Template.txt \
		-C -R --ALIGNERS "gmap,blat" -N 1 \
		--MAX_INTRON_LENGTH ${max_intron_len} ${cuff} --CPU 20 \
		-g pasa/genome.fa -t pasa/transcripts.fasta --TDN pasa/trinity_denovo.accs \
		${orient} -L --annots_gff3 ${gff3} \
		--gene_overlap 50
fi

echo "Generate the comprehensive transcriptome database -- `date`" >&2
${pasa_home}/scripts/build_comprehensive_transcriptome.dbi \
	-c pasa/pasa.alignAssembly.Template.txt -t pasa/transcripts.fasta \
	--min_per_ID 95 --min_per_aligned 30

echo "Finished $0 $@ -- `date`" >&2


