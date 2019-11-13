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

This script is used to upload files to an FTP server using lftp.

${txtbld}OPTIONS${txtrst}:
	-S	Start [Giving any string to start upload]
	-f	FTP address ${bldred}[NECESSARY,  default ftp://submit.big.ac.cn]${txtrst}
	-u	User name ${bldred}[NECESSARY,  default chentong_biology@163.com]${txtrst}
	-p	Password ${bldred}[NECESSARY, default 550336392]${txtrst}
	-t	Target dir ${bldred}[NECESSARY,  default GSA]${txtrst}
	-s	Source dir ${bldred}[NECESSARY, default upload]${txtrst}	
	-T	BIG database
EOF
}

start_now=
ftp=ftp://submit.big.ac.cn
user=chentong_biology@163.com
passwd=550336392
target=GSA
source_dir=upload

while getopts "hf:u:p:t:S:s:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		S)
			start_now=$OPTARG
			;;
		f)
			ftp=$OPTARG
			;;
		u)
			user=$OPTARG
			;;
		p)
			passwd=$OPTARG
			;;
		t)
			target=$OPTARG
			;;
		s)
			source_dir=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z $start_now ]; then
	usage
	exit 1
fi

cat <<END >lftp.script
open -u ${user},${passwd} ${ftp}
#mkdir -p ${target}
cd ${target}
cache size 33554432
set cmd:parallel 10
mput -c ${source_dir}/*
END

lftp -f lftp.script
