#!/bin/bash

#set -x

filename=`basename $0`

usage()
{
cat <<EOF
${txtcyn}
Usage:

${filename} options${txtrst}

${bldblu}Function${txtrst}:

This is a collection of FASTA tools.

#### General
${filename} computeStatPandaMatrix
${filename} summaryMatrix
${filename} summaryMatrix2
${filename} normPandaMatrix
${filename} pasteMultipleFilesPandaCustom
${filename} pasteMultipleFilesPanda
${filename} pearsonCorrelationMatrix
${filename} pearsonCorrelation
${filename} removeSepcialColumnsFromMatrix
${filename} extractSpecialColumnsFromMatrix
${filename} methylC_concat
${filename} pasteMethylC
${filename} transferPandasSummaryForBoxplotUsage
${filename} renameColumnsFromMatrix
${filename} extractSpecialRowsFromMatrix
${filename} filterStatPandaMatrixByGroup

#### Special
${filename} combineEncodeData
${filename} combineTCGAData
${filename} concatEncodeData (faster)
${filename} concatTCGAData (faster)
${filename} concatTCGAData2

EOF
}

if test $# -lt 1; then
	usage
	exit 1
fi

program=$1".py"
type ${program} >/dev/null 2>&1
error=$?
if test $error != 0; then
	usage
	echo "**Please check the program name input**"
	exit 1
else
	shift
	${program} "$@"
fi

