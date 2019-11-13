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

This is designed to run generate program for special environment uages. 

Using this script, one can just like using normal command without realizing the inner environment activation, program run and environment deactivation.

${txtbld}OPTIONS${txtrst}:
	-c	Program name without path ${bldred}[NECESSARY]${txtrst}
	-e	Environment name${bldred}[NECESSARY]${txtrst}
	-b	Conda path${bldred}[Optional. Default: /miniconda2/bin]${txtrst}
	-o	Output file name. Normally [/MPATHB/bin/<Program name>]
EOF
}

command_cmd=''
environment=''
conda_path='/miniconda2/bin'
output_cmd=''

while getopts "hc:e:b:o:" OPTION
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
		o)
			output_cmd=$OPTARG
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

#if ! [ -z ${conda_path} ]; then
#	export PATH=${conda_path}:${PATH}
#fi

if test -s ${output_cmd}; then
	(>&2 echo "***${output_cmd} exists, Please check first.***")
	exit 1
fi

cat <<END >${output_cmd}
#!/bin/bash
# set -x # for debugging
source ${conda_path}/activate ${environment}
${command_cmd} \$*
source ${conda_path}/deactivate ${environment}
END

chmod 755 ${output_cmd}
