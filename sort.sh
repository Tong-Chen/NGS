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

This script is used to do perform sort operation with header line kept.

${txtbld}OPTIONS${txtrst}:
	-f	Data file ${bldred}[NECESSARY, - represents STDIN]${txtrst}
	-p	Parameters for sort command${bldred}[NECESSARY, lowercase p]${txtrst}
	-z	Number of header lines[${bldred}Default 1${txtrst}]
EOF
}

file=
sortP=
header=1

while getopts "hf:p:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		p)
			sortP=$OPTARG
			;;
		z)
			header=$OPTARG
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

if test "$file" == "-"; then
	cat $file >___tmp.59
	file="___tmp.59"
fi

head -n ${header} ${file:-/dev/stdin}

tail -n +$(( ${header}+1 )) ${file:-/dev/stdin} | sort ${sortP}

/bin/rm -f ___tmp.59
