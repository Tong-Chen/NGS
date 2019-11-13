#!/bin/bash

echo "This lists the information of this computer."

echo

echo "Hostname is $(tput setaf 3)`hostname`$(tput sgr0),\
Ip address is $(tput setaf 3)\
`/sbin/ifconfig | sed -n '2p' | cut -d ':' -f 2 | cut -d ' ' -f 1`.
$(tput sgr0)"


nuclear=`uname -a | cut -d ' ' -f 3`
bitInfo=`uname -a | cut -d ' ' -f 12`

if test $bitInfo == "x86_64"; then
    bit=64
else
    bit=32
fi

echo "The $(tput bold)${bit}$(tput sgr0) bit operating \
system is $(tput bold) `head -n 1 /etc/issue`\
$(tput sgr0), Nuclear info is $(tput setaf 1)\
${nuclear}$(tput sgr0)."

echo

echo "The CPU is$(tput setaf 4)`sed -n '5p' /proc/cpuinfo \
| cut -d ':' -f 2 | sed 's/[ ] */ /g'`$(tput sgr0)."

echo

echo "There are $(tput setaf 5)\
`cat /proc/cpuinfo | grep "physical id" | sort | uniq \
| wc -l`$(tput sgr0) physical cpu, \
each physical \
cpu has$(tput setaf 5)`sed -n '12p' /proc/cpuinfo | \
cut -d ':' -f 2`$(tput sgr0) cores,\
$(tput setaf 5)`sed -n '10p' /proc/cpuinfo | \
cut -d ':' -f 2`$(tput sgr0) threads."

echo

echo "There are $(tput setaf 5)\
`cat /proc/cpuinfo | grep "cpu cores" | wc -l`$(tput sgr0) logical cpu."

mem=`head -n 1 /proc/meminfo | cut -d ':' -f 2 | sed 's/^ *//g' | cut -d ' ' -f 1`
memInM=$(echo "$mem/1024/1024" | bc -l)

echo

echo "The memory of this server is $(tput setaf 5)${memInM}$(tput sgr0)G."

echo

echo "The disk information is :"

echo "`df -h`"

