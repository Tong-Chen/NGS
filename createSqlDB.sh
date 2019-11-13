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

This script is used to create mysql databases.

${txtbld}OPTIONS${txtrst}:
	-f	Database name ${bldred}[NECESSARY]${txtrst}
	-r	Root name or other username with privideleges to create
		database.
	#-p	Password for username given to <-r>.
	-u	Username for operations on new created database.
	-q	Password for username given to <-u>. 
		Password must contain number, alphabets in both lowercase and uppercase, 
	    special symbols.	
EOF
}

dbname=
header='TRUE'

while getopts "hf:r:p:u:q:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			dbname=$OPTARG
			;;
		r)
			root_name=$OPTARG
			;;
		p)
			root_passwd=$OPTARG
			;;
		u)
			user_name=$OPTARG
			;;
		q)
			user_passwd=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z $dbname ]; then
	usage
	exit 1
fi

cat <<END | mysql -u ${root_name} -p
CREATE DATABASE ${dbname};
GRANT SELECT,INSERT,UPDATE,DELETE,CREATE VIEW,CREATE,INDEX,DROP on ${dbname}.* TO '${user_name}'@'localhost' IDENTIFIED BY '${user_passwd}';
FLUSH PRIVILEGES;
END
