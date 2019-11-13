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

This script is used to do ********************.

${txtbld}OPTIONS${txtrst}:
	-u	Username ${bldred}[NECESSARY]${txtrst}
	-f	First-time create ${bldred}[FALSE]${txtrst}
EOF
}

username=
first_time='FALSE'

while getopts "hu:f:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		u)
			username=$OPTARG
			;;
		f)
			first_time=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z $username ]; then
	usage
	exit 1
fi

if [ ${first_time} == "TRUE" ]; then
	groupadd sftpusers
	cp /etc/ssh/sshd_config /etc/ssh/sshd_config.bak
	sed -i 's/^Subsystem/#Subsystem/' /etc/ssh/sshd_config
	cat <<END >>/etc/ssh/sshd_config
Subsystem       sftp    internal-sftp
Match Group sftpusers
ChrootDirectory /MPATHD/ftp/%u
ForceCommand    internal-sftp
AllowTcpForwarding no 
X11Forwarding no
END
	cp /etc/selinux/config /etc/selinux/config.bak
	sed -i 's/SELINUX=enforcing/SELINUX=disabled/' /etc/selinux/config
	setenforce 0 
fi

useradd -s /bin/false -G sftpusers ${username}
passwd ${username}

mkdir -p /MPATHD/ftp/${username}/upload
usermod -d /MPATHD/ftp/${username} ${username}


chown root:sftpusers /MPATHD/ftp/${username}
chmod 755 /MPATHD/ftp/${username}

chown ${username}:sftpusers /MPATHD/ftp/${username}/upload
chmod 755 /MPATHD/ftp/${username}/upload

service sshd restart



