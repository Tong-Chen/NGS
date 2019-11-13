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

This script is used to do ********************.

${txtbld}OPTIONS${txtrst}:
	-f	File names without <.pdb> extension ${bldred}[NECESSARY]${txtrst}
		Multiple names separated by one blank are accepted.
	-t	Analysisi type <procheck>(default) or <procheck_cmp>.
	-o	Output dir for compare results.
	-r	The resolution at which the structure has been determined. 
		If the resolution is unknown or not relevant 
		(eg for structures solved by NMR) enter any resolution (eg 2.0). 
		${bldred}[NECESSARY]${txtrst}
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
EOF
}

file=
header='TRUE'
resolution=2.0
type_l='procheck'
out=

while getopts "hf:r:t:o:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		t)
			type_l=$OPTARG
			;;
		o)
			out=$OPTARG
			;;
		r)
			resolution=$OPTARG
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

if [ -z "$file" ]; then
	usage
	exit 1
fi

if [ "${type_l}" == "procheck" ]; then
	for eachfile in ${file}; do
		mkdir -p ${eachfile}.procheck
		/bin/cp -f ${eachfile}.pdb ${eachfile}.procheck
		(cd ${eachfile}.procheck; procheck ${eachfile}.pdb ${resolution})
		(cd ${eachfile}.procheck; for i in ${eachfile}_*.ps; do ps2pdf $i ${i/ps/pdf}; done)
		(cd ${eachfile}.procheck; /bin/rm -f *.ps)
	done
fi


if [ "${type_l}" == "procheck_comp" ]; then
	mkdir -p ${out}
	echo "${file}" | tr ' ' '\n' >${out}/${out}
	(cd ${out}; procheck_comp ${out}/${out})
	(cd ${out}; for i in ${out}_*.ps; do ps2pdf $i ${i/ps/pdf}; done)
	(cd ${out}; /bin/rm -f *.ps)
fi
