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
  
This script is used to upload files to an FTP server using lftp.
  
${txtbld}OPTIONS${txtrst}:
	-f	FTP address ${bldred}[NECESSARY]${txtrst}
	-u	User name ${bldred}[NECESSARY]${txtrst}
	-p	Password ${bldred}[NECESSARY]${txtrst}
	-t	Target dir ${bldred}[NECESSARY, for GEO in format like
		<fasp/GEO_metadata_zhaohui>]${txtrst}
	-s	Source dir ${bldred}[NECESSARY, default current directory]${txtrst}	
	-r	Send success information to given email.
		${bldred}[OPTIONAL, default chentong_biology@163.com]${txtrst}	
EOF
}
  
ftp=
user=
passwd=
target=
source_dir="."
email="chentong_biology@163.com"
  
while getopts "hf:u:p:t:s:r:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
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
		r)
			email=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done
  
if [ -z $ftp ]; then
	usage
	exit 1
fi
  
cat <<END >lftp.script
open -u ${user},${passwd} ${ftp}
mkdir -p ${target}
cd ${target}
cache size 33554432
set cmd:parallel 10
mput -c ${source_dir}/*
END
  
lftp -f lftp.script

pythonmail3.py -r ${email} -s "${target} FINISHED"

