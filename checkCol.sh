#!/bin/bash

#set -x
set -e 
set -u

usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to check the number of columns in given line and
will output the column number before each column item.

${txtbld}OPTIONS${txtrst}:
	-f	Data file (- represents STDIN)${bldred}[NECESSARY]${txtrst}
	-l	lineNumber[${txtred}Default 1 meaning the first line${txtrst}]
EOF
}

file=
lineno=1

while getopts "hf:l:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		l)
			lineno=$OPTARG
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

if read -t 0; then
	sed -n "${lineno} s/\t/\n/g; ${lineno}p" | \
		sed '=' | sed 'N;s/\n/\t/'
else
	sed -n "${lineno} s/\t/\n/g; ${lineno}p" ${file} | \
		sed '=' | sed 'N;s/\n/\t/'
fi
