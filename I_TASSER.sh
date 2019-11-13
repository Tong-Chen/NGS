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

This script is used to facilitate I-TASSER usage.
Output will be directory given to <-d>.

${txtbld}OPTIONS${txtrst}:
	-f	Directory for I-TASSER${bldred}[NECESSARY]${txtrst}
		Default: /MPATHB/soft/I_TASSER/I-TASSER5.0/
	-l	Directory for template libraries${bldred}[NECESSARY]${txtrst}
		Default: /MPATHB/soft/I_TASSER/ITLIB/
	-j	Directory containing "bin/java"${bldred}[NECESSARY]${txtrst}
		Default: /usr/
	-n	Name of your query sequence.${bldred}[NECESSARY]${txtrst}
	-d	Directory containing file <seq.fasta> which saves query sequence 
		(without terminating * at the end of amino acid sequences)
	-o	Specify one or all other analyses like
		<-LBS>: predict ligand-binding site[default]
		<-EC>:	predict EC number
		<-GO>:	predict GO terms
	-H	Specify the maximum hours of simulations. 
		${bldred}[Default 20]${txtrst}
	-F	Specify to run quick simulations.
		${bldred}[default TRUE,  which used <-light parameter>]${txtrst}
EOF
}

pkgdir=/MPATHB/soft/I_TASSER/I-TASSER5.0
libdir=/MPATHB/soft/I_TASSER/ITLIB/
javadir=/usr/
seq_name=
seq_dir=
other_par="-LBS"
hours=20
light=TRUE

while getopts "hf:l:j:n:d:o:H:F" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			pkgdir=$OPTARG
			;;
		l)
			libdir=$OPTARG
			;;
		j)
			javadir=$OPTARG
			;;
		n)
			seq_name=$OPTARG
			;;
		d)
			seq_dir=$OPTARG
			;;
		o)
			other_par=$OPTARG
			;;
		H)
			hours=$OPTARG
			;;
		F)
			light=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z ${seq_name} ]; then
	usage
	exit 1
fi

if test "${light}" == "TRUE"; then
	light="-light"
else
	light=""
fi
cat <<EOF
perl ${pkgdir}/I-TASSERmod/runI-TASSER.pl -pkgdir ${pkgdir} -libdir ${libdir} \
	-java ${javadir} -runstyle gnuparallel \
	-seqname ${seq_name} -datadir ${seq_dir} \
	${other_par} ${light} -hours ${hours}
EOF


perl ${pkgdir}/I-TASSERmod/runI-TASSER.pl -pkgdir ${pkgdir} -libdir ${libdir} \
	-java ${javadir} -runstyle gnuparallel \
	-seqname ${seq_name} -datadir ${seq_dir} \
	${other_par} ${light} -hours ${hours}

	#-homoflag benchmark 
