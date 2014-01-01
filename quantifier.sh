#!/bin/bash

set -x

usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to quantify miRNA exoression using `quantifier.pl`

${txtbld}OPTIONS${txtrst}:
	-p	precursor.fa  miRNA precursor sequences from miRBase
 	-m	mature.fa     miRNA sequences from miRBase
	-r	reads.fa      your read sequences
	-y	output_suffix 
	-o	Other parameters given to a modified 
		version of quantifier.pl <quantifier.modified.pl>.
		Default '-j -W'.
	-O	Additional parameters other than those 
   		given to <-o> given to a modified 
		version of quantifier.pl <quantifier.modified.pl> 
		in triming stage. 
		Default '-U -C 20'.
		All parameters given to -o and -O will be used in trimming
		stage,  but only those given to -o will be sued at first round
		of quantification.
	-t	[species] e.g. Mouse or mmu], if not searching in a specific
		species all species in your files will be analyzed else only the
		species in your dataset is considered.
	-T	Trim unmapped reads	[Default TRUE]
	-P	Parameters given to <trimFasta.py>. 
		Default "-s 0 -e -1 -l 16"
	-n	Number of trimming you want to perform.
		Default perfom trimming until 
		empty files generated. Accept a number to indicate the
		maximum allowed trim time.
EOF
}

precursor=
mature=''
reads=''
suffix=''
spe=''
other_p_quantifier='-j -W'
additional_p_quantifier="-U -C 20"
trim="TRUE"
parameter_trimFasta="-s 0 -e -1 -l 16"
num_trim='A'

while getopts "h:p:m:r:y:o:t:T:P:n:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		p)
			precursor=$OPTARG
			;;
		m)
			mature=$OPTARG
			;;
		r)
			reads=$OPTARG
			;;
		y)
			suffix=$OPTARG
			;;
		t)
			spe=$OPTARG
			;;
		o)
			other_p_quantifier=$OPTARG
			;;
		O)
			additional_p_quantifier=$OPTARG
			;;
		T)
			trim=$OPTARG
			;;
		P)
			parameter_trimFasta=$OPTARG
			;;
		n)
			num_trim=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done
if [ -z $precursor ] || [ -z $mature ] || [ -z $reads ] || [ -z $suffix ]; then
	usage
	exit 1
fi

quantifier.modified.pl -p $precursor -m ${mature} \
	-r ${reads} -y ${suffix} -t ${spe} ${other_p_quantifier}

if test "${trim}" == 'TRUE'; then
	#---Get unmapped reads
	getSeqOfLocus.py -i ${reads} -l \
	    <(cut -f 1 expression_analyses/expression_analyses_${suffix}/read_occ) \
		-o exclude >${reads}.unmapped.fa
	#---Get mapped reads
	getSeqOfLocus.py -i ${reads} -l \
	    <(cut -f 1 expression_analyses/expression_analyses_${suffix}/read_occ) \
		-o include >${reads}.mapped.fa
	#---------Trim unmapped reads------------------
	trimFasta.py -i ${reads}.unmapped.fa ${parameter_trimFasta} \
		>${reads}.unmapped.trimed.fa
	trimTime=0
	while test -s ${reads}.unmapped.trimed.fa; do
		#---New reads-------------
		#cat ${reads}.mapped.fa ${reads}.unmapped.trimed.fa \
		#	>${reads}.${suffix}.processed.fa
		cp -f ${reads}.unmapped.trimed.fa ${reads}.processed.fa
		/bin/rm -rf expression_analyses/expression_analyses_${suffix} \
			miRNAs_expressed_all_samples_${suffix}
		#---Test unmapped trimmed reads----------
		quantifier.modified.pl -p $precursor -m ${mature} \
			-r ${reads}.processed.fa -y ${suffix} -t ${spe} \
			${other_p_quantifier} \
			${additional_p_quantifier}
		#---Get unmapped reads
		getSeqOfLocus.py -i ${reads}.processed.fa -l \
			<(cut -f 1 expression_analyses/expression_analyses_${suffix}/read_occ) \
			-o exclude >${reads}.unmapped.fa
		#---Get mapped reads
		getSeqOfLocus.py -i ${reads}.processed.fa -l \
			<(cut -f 1 expression_analyses/expression_analyses_${suffix}/read_occ) \
			-o include >>${reads}.mapped.fa
		#---------Trim unmapped reads------------------
		((trimTime++))
		echo "--Perfomed "${trimTime}" trims--"
		if test ${num_trim} != "A" && test ${trimTime} -ge ${num_trim}; then 
			break
		fi
		trimFasta.py -i ${reads}.unmapped.fa ${parameter_trimFasta} \
			>${reads}.unmapped.trimed.fa
	done
	echo "--Perform final quantification---"
	/bin/rm -rf expression_analyses/expression_analyses_${suffix} \
		miRNAs_expressed_all_samples_${suffix}
	quantifier.modified.pl -p $precursor -m ${mature} \
		-r ${reads}.mapped.fa -y ${suffix} \
		-t ${spe} ${other_p_quantifier}
fi
