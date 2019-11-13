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

This script is used to generate documents.

${txtbld}OPTIONS${txtrst}:
	-t	Title for each section
		${bldred}[NECESSARY, like "Sequencing quality"]
		The supplied text will be treated as h1 and will start a new
		section.
		${txtrst}
	-l	List of files to be cat together.
EOF
}

title=
file_list=

while getopts "ht:l:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		t)
			title=$OPTARG
			;;
		l)
			file_list=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z "$title" ]; then
	usage
	exit 1
fi

cat <<END
# ${title}

END

if [ -n "${file_list}" ]; then
	cat ${file_list}
fi
