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

This script is used to perform orthoMcl analysis using MySql, MCL and
orthomcl.

Before running this script, one must have one mysql database and a
mysql user which can perform operation on this database.


${txtbld}OPTIONS${txtrst}:
	-d	Mysql database name (using user_name as prefix to avoid
		duplication) ${bldred}[Necessary]${txtrst}
	-u	Mysql database username ${bldred}[Necessary]${txtrst}
	-p	Mysql database password ${bldred}[Necessary]${txtrst}
	-s	Target species of this analysis 
		(Any representing string is OK, the shorter the better)
		${bldred}[Necessary]${txtrst}
	-D	A directory containing FASTA files for all proteins.
		${bldred}[Necessary]${txtrst}
	-S	Sequences downloaded from orthMCL website.
		${bldred}[Optional, for example:
		/MPATHB/data/orthomcl/aa_seqs_OrthoMCL-5.fasta]${txtrst}
	-t	Number of threads for blast. ${bldred}[Default 50]${txtrst}
EOF
}

database=
user=
passwd=
prefix=
datadir=
orthofa=
threads=50

while getopts "hd:u:p:s:D:S:t:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		d)
			database=$OPTARG
			;;
		u)
			user=$OPTARG
			;;
		p)
			passwd=$OPTARG
			;;
		s)
			prefix=$OPTARG
			;;
		D)
			datadir=$OPTARG
			;;
		S)
			orthofa=$OPTARG
			;;
		t)
			threads=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z $database ]; then
	usage
	exit 1
fi

cat <<END >orthomcl.config
dbVendor=mysql 
dbConnectString=dbi:mysql:${database}
dbLogin=${user}
dbPassword=${passwd}
# Change strings as you like
similarSequencesTable=${prefix}SimilarSequences
orthologTable=${prefix}Ortholog
inParalogTable=${prefix}InParalog
coOrthologTable=${prefix}CoOrtholog
#Standards
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE
END

if ! test -s ${prefix}.orthoInit.OK; then
	orthomclInstallSchema orthomcl.config inst_schema.log
	echo "${prefix}.orthoInit.OK" >${prefix}.orthoInit.OK
fi

if ! test -s ${prefix}.goodProt.fasta; then
	orthomclFilterFasta ${datadir} 10 20 \
		${prefix}.goodProt.fasta ${prefix}.poorProt.fasta 
fi

if ! test -s ${prefix}.orthl.fasta; then
	cat ${prefix}.goodProt.fasta ${orthofa} >${prefix}.orthl.fasta
fi
	
#makeblastdb -in ${prefix}.orthl.fasta -dbtype prot -title orthomcl \
#	-parse_seqids -out ${prefix}.orthomcl -logfile orthomcl.log

if ! test -s ${prefix}.orthomcl.psq; then
	makeblastdb -in ${prefix}.orthl.fasta -dbtype prot -title orthomcl \
		-out ${prefix}.orthomcl -logfile orthomcl.log
fi

if ! test -s ${prefix}.orthomcl.blastout; then
	blastp -db ${prefix}.orthomcl -query ${prefix}.goodProt.fasta -seg yes \
		-out ${prefix}.orthomcl.blastout -evalue 1e-5 -outfmt 7 \
		-num_threads ${threads}
fi

if ! test -s ${prefix}.orthomcl.blastout.2; then
	grep -v '^#' ${prefix}.orthomcl.blastout >${prefix}.orthomcl.blastout.2
fi

if ! test -s ${prefix}.similarSequences.txt; then
	orthomclBlastParser ${prefix}.orthomcl.blastout.2 ${datadir} \
		>${prefix}.similarSequences.txt
	perl -p -i -e 's/0\t0/1\t-181/' ${prefix}.similarSequences.txt
fi

#perl -p -i -e 's/\t(\w+)(\|.*)orthomcl/\t$1$2$1/' ${prefix}.similarSequences.txt

if ! test -s ${prefix}.orthomclLoadBlast.OK; then
	orthomclLoadBlast orthomcl.config ${prefix}.similarSequences.txt
	echo "${prefix}.orthomclLoadBlast.OK" >${prefix}.orthomclLoadBlast.OK
fi

if ! test -s ${prefix}.orthomclPairs.OK; then
	orthomclPairs orthomcl.config orthomcl_pairs.log cleanup=no
	echo "${prefix}.orthomclPairs.OK" >${prefix}.orthomclPairs.OK
fi

if ! test -s mclInput; then
	orthomclDumpPairsFiles orthomcl.config
fi

if ! test -s mclOutput; then
	mcl mclInput --abc -I 1.5 -o mclOutput -te ${threads}
fi

if ! test -s ${prefix}.groups.xls; then
	orthomclMclToGroups ${prefix}_ 10000 <mclOutput >${prefix}.groups.xls
fi

