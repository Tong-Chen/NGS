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
	-u	User name ${bldred}[NECESSARY]${txtrst}
	-p	Password ${bldred}[NECESSARY]${txtrst}
	-g	Group ${bldred}[Currently support ehbio, probation, intern]${txtrst}
EOF
}

user=
password=
group=

while getopts "hu:p:g:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		u)
			user=$OPTARG
			;;
		p)
			password=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z $user ]; then
	usage
	exit 1
fi

cat <<END >/tmp/${user}.create
${user}:${password}::${group}::/MPATHC/home/${user}:/bin/bash
END

newusers </tmp/${user}.create


/bin/cp -u /MPATHC/home/.bashrc /MPATHC/home/${user}/.bashrc
chown ${user}:${group} /MPATHC/home/${user}/.bashrc
/bin/cp -u /MPATHC/home/.bash_profile /MPATHC/home/${user}/.bash_profile
chown ${user}:${group} /MPATHC/home/${user}/.bash_profile
chmod go-rwx /MPATHC/home/${user}


