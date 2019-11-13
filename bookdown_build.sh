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

This script is used to build bookdown html or PDF.

index.Rmd is needed.

${txtbld}OPTIONS${txtrst}:
	-w	Build HTML ${bldred}[Default FALSE]${txtrst}
		Accept a file name like <index.Rmd> or <index_html.Rmd> to start build web from this file.
	-p	BUild PDF ${bldred}[Default FALSE]${txtrst}
		Accept a file name like <index.Rmd> or <index_pdf.Rmd> to start build PDF from this file.
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
EOF
}

html='FALSE'
pdf='FALSE'
header='TRUE'

while getopts "hw:p:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		w)
			html=$OPTARG
			;;
		p)
			pdf=$OPTARG
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

if [ $html == "FALSE" ] && [ $pdf == "FALSE" ]; then
	usage
	exit 1
fi

if [ $html != "FALSE" ]; then
	Rscript -e "bookdown::render_book('$html','bookdown::gitbook')"
fi

if [ $pdf != "FALSE" ]; then
	Rscript -e "bookdown::render_book('$pdf','bookdown::pdf_book')"
fi


