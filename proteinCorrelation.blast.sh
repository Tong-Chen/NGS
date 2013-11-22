#!/bin/bash

set -e 
set -u

usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to do protein correlation analysis and protein
protein interaction sites using clustalw, protdist(phylip) and DSSP.

${txtbld}OPTIONS${txtrst}:
	-f	Data file1 (Fasta, protein 1)
		${bldred}[NECESSARY, y axis]${txtrst}
	-g	Data file2 (Fasta, protein 2)
		${bldred}[NECESSAR, x axis]${txtrst}
	-b	Blast filter identity.
		[Default: larger than 0.7, accept a float number]
	-d	Length deviation.
		${bldred}[Default less than 0.2, accept a float number]${txtrst}
	-c	Evalue for blast output.[default blastp default 10]
	-a	Delete formula.
		${bldred}[Default 2, accept 0 or 1. 
		0 means not delete same sequences. 1 means delete sequences
		same in one dataset. 2 means delete sequences same in both
		datasets.]${txtrst}
	-e	Species of aim sequence.
		${bldred}[NECESSARY when -d is a number, use for getting protein 
		length for ]${txtrst}
	-t	Threshold for correlation[${txtred}Necessary, default 0${txtrst}]
	-p	Compute PPI sites[${txtred}Default FALSE${txtrst}]
	-s	The window size for computing PPI sites.
		[${txtred}Default 5 whne -p is TRUE${txtrst}]
	-u	The sliding distance of each window.
		[${txtred}Default same with -s whne -p is TRUE${txtrst}]
	-q	Output contour picture.[${txtred}Default FALSE${txtrst}]

General usage:

proteinCorrelation.blast.sh -f gp5.fa -g gp5.2.fa -p TRUE -q TRUE
EOF
}

file1=
file2=
del_formula=2
identity=0.7
lendev=0.2
species=
thresh=0
ppi='FALSE'
size=5
slide=$size
contour='FALSE'
nr='/data7T/mercu-b-backup/pub-data/NCBI-NR/nr'
taxid='~/home/server-project/ncbi/gi_taxid_prot.dmp'
evalue=10

while getopts "hf:g:a:b:d:c:e:t:p:s:q:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file1=$OPTARG
			;;
		g)
			file2=$OPTARG
			;;
		a)
			del_formula=$OPTARG
			;;
		b)
			identity=$OPTARG
			;;
		d)
			lendev=$OPTARG
			;;
		c)
			evalue=$OPTARG
			;;
		e)
			species=$OPTARG
			;;
		t)
			thresh=$OPTARG
			;;
		p)
			ppi=$OPTARG
			;;
		s)
			size=$OPTARG
			;;
		u)
			slide=$OPTARG
			;;
		q)
			contour=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done
if [ -z $file1 ] || [ -z $file2 ]; then
	usage
	exit 1
fi

if [ "$file1" = "$file2" ]; then
	echo "Two files with same name are not allowed!"
	exit 1
fi

mid=${size}.${slide}

phylip_dist ()
{
	#ln -sf $1 infile
	#/bin/rm -f outfile
	echo 'Y' | protdist $1 $1.Dist 
	#>/dev/null
	#/bin/mv -f outfile $1.Dist
	#/bin/rm -f infile
}

align_dist ()
{
	clustalw $1 -TYPE=PROTEIN -OUTPUT=CLUSTAL \
		-OUTFILE=$1.aln -OUTORDER=INPUT >$1.log
	clustalw2phylip $1.aln >$1.Phylip
	phylip_dist $1.Phylip >/dev/null
	echo ">>>Finish $1"
}

split_align_dist()
{
	/bin/rm -f `ls | grep -P "^$1\.[0-9]+$"`
	splitPhylip.py $1 $size $slide
	for i in `ls | grep -P "^$1\.[0-9]+$"`; do \
		phylip_dist $i >/dev/null & done
	wait
}

blast()
{
	blastp -query $1 -db ${nr} -num_threads=8 -evalue ${evalue} \
		-outfmt='6 qseqid sseqid pident length slen evalue' \
		| awk -v iden="${identity}" -v ld="${lendev}" \
		'BEGIN{OFS="\t";FS="\t"}{diff=$5-$4; \
		diff=(diff>0)?diff:(-1)*diff; \
		if( $3>=iden && diff/$3<=ld) print $0}' \
		| tee $1_out \
		| cut -f 2 | cut -d '|' -f 2 > $1_gi
	awk 'BEGIN{OFS="\t";FS="\t"}{if($3==100.0 && $4==$5) print $2}' \
		$1_out | cut -f 2 -d '|' -s >$1_aim.gi &
	cat <<END >$1_sql
drop table if exists $2;
create table $2 (gi INT, PRIMARY KEY (gi));
load data local infile "`pwd`/$1_gi" into table $2;
select gi_taxid_prot.gi, gi_taxid_prot.taxid from gi_taxid_prot,$2 where $2.gi = gi_taxid_prot.gi; 
drop table if exists $2;
END
	mysql -h 210.75.224.29 -upc -ppc ncbi_taxonomy \
		<$1_sql >$1_taxid &
	blastdbcmd -entry_batch $1_gi -db ${nr} -dbtype=prot \
		>$1_hit &
	wait
	if test -e $1_aim.gi && (! test -s $1_aim.gi); then
		echo "****Warning: $1 is not found in NCBI database.****"
	fi	
	diffComByOneCol.py $1_aim.gi $1_taxid >$1_aim.gi_taxid
	parseBlastOutputForPC.py $1_out $1_hit $1_taxid >$1_set

	echo ">>>Finish $1"
}

blast $file1 tablea &
blast $file2 tableb &

oldfile1=$file1
file1=${file1}_set
oldfile2=$file2
file2=${file2}_set

wait

echo ">>>Finish blast"

preprocessFastaForPorteinCorrelation.py $file1 $file2 ${del_formula}

if [ $(cat ${file1}_${file2}.Rev | wc -l) -lt 4 ] || \
   [ $(cat ${file2}_${file1}.Rev | wc -l) -lt 4 ]
then
	echo ">>>Less sequences."
	exit 1
fi


align_dist ${file1}_${file2}.Rev &
align_dist ${file2}_${file1}.Rev &

wait

phylipDistanceMatrix.py ${file1}_${file2}.Rev.Phylip.Dist \
	${file2}_${file1}.Rev.Phylip.Dist \
	`cat ${file1}_${file2}.Rev.Phylip.Dist | wc -l` | tee \
	>${file1}_${file2}.Correlation
	
if test "$ppi" == "TRUE"; then
	(split_align_dist ${file1}_${file2}.Rev.Phylip &
	split_align_dist ${file2}_${file1}.Rev.Phylip &
	wait)
	/bin/rm -f ${file1}_${file2}.Correlation.Matrix
	for i in `/bin/ls | /bin/grep -P "^${file1}_${file2}.Rev.Phylip.[0-9]+.Dist$"`; do
		label1=${i/${file1}_${file2}.Rev.Phylip./}
		label1=${label1/.Dist/}
		no=0
		for j in `/bin/ls | /bin/grep -P "^${file2}_${file1}.Rev.Phylip.[0-9]+.Dist$"`; do 
			no=$(( no+1 ))
			label2=${j/${file2}_${file1}.Rev.Phylip./}
			label2=${label2/.Dist/}
			/bin/echo -e \
				"$label1""\t""$label2""\t"$(phylipDistanceMatrix.py $i $j `cat $i | wc -l`) \
				>${file1}_${file2}.Correlation.Matrix.${no} &
				#>>${file1}_${file2}.Correlation.Matrix
		done
		wait
		cat `ls | grep -P \
			"^${file1}_${file2}\.Correlation\.Matrix\.[0-9]+$"` \
			>>${file1}_${file2}.Correlation.Matrix
		/bin/rm -f `ls | grep -P \
			"^${file1}_${file2}\.Correlation\.Matrix\.[0-9]+$"`
	done

fi


if test "$ppi" == 'TRUE' && test "$contour" == 'TRUE'; then 
	transformCorMatrix.py ${file1}_${file2}.Correlation.Matrix \
		>${file1}_${file2}.Correlation.Matrix.Plot
	heatmap.2.sh -f ${file1}_${file2}.Correlation.Matrix.Plot \
		-k FALSE -r FALSE -b '0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1'
fi

#------------mv things--------------------------------
/bin/rm -rf ${oldfile1}.${oldfile2}.${mid}
/bin/mkdir -p ${oldfile1}.${oldfile2}.${mid}
#/bin/mv ${file1}.* ${file1}_${file2}
#/bin/mv ${file2}.* ${file1}_${file2}
/bin/mv ${oldfile1}_* ${oldfile1}.${oldfile2}.${mid}
/bin/mv ${oldfile2}_* ${oldfile1}.${oldfile2}.${mid}
#/bin/mv ${file1}_* ${file1}.${file2}
#/bin/mv ${file2}_* ${file1}.${file2}
/bin/cp -f ${oldfile1}.${oldfile2}.${mid}/*.png .
/bin/cp -f ${oldfile1}.${oldfile2}.${mid}/*.aln .
