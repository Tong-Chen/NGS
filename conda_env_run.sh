#!/bin/bash

#set -x

usage()
{
cat <<EOF
${txtcyn}

***CREATED BY Chen Tong (chentong_biology@163.com)***

Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This is designed to run conda program in given environment. It will automatically activate the environment, run the program and deactivate the environment.

Thress commands from python, 'activate', 'conda', 'deactivate' must be in PATH.

${txtbld}OPTIONS${txtrst}:
	-c	Full command to be run ${bldred}[NECESSARY]${txtrst}
	-e	Environment name${bldred}[NECESSARY]${txtrst}
	-b	Conda path${bldred}[Optional. Default: /miniconda2/bin]${txtrst}
EOF
}

command_cmd=''
environment=''
conda_path='/miniconda2/bin'

while getopts "hc:e:b:" OPTION
do
	case $OPTION in
		h)
			echo "Help mesage"
			usage
			exit 1
			;;
		c)
			command_cmd=$OPTARG
			;;
		e)
			environment=$OPTARG
			;;
		b)
			conda_path=$OPTARG
			;;
		?)
			usage
			echo "Unknown parameters"
			exit 1
			;;
	esac
done


if [ -z ${environment} ]; then
	echo 1>&2 "Please give command and environment."
	usage
	exit 1
fi

if ! [ -z ${conda_path} ]; then
	export PATH=${conda_path}:${PATH}
fi

source activate ${environment}
${command_cmd}
source deactivate ${environment}

