#!/bin/bash

#set -x
set -e
#set -u

usage()
{
cat <<EOF >&2
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to transfer GTF to bed12 format.

${txtbld}OPTIONS${txtrst}:
	-f	GTF file ${bldred}[NECESSARY]${txtrst}
EOF
}

file=
header='TRUE'

while getopts "hf:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
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

if [ -z $file ]; then
	usage
	exit 1
fi


gtfToGenePred -ignoreGroupsWithoutExons ${file} GRCh38.gtf.50505050.pred
genePredToBed GRCh38.gtf.50505050.pred ${file}.bed12
/bin/rm -f GRCh38.gtf.50505050.pred
