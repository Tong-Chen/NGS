#!/bin/bash


usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to correct wrong tags from Tophat output.
Newly corrected file has the format *.corrected.bam ot *.corrected.sam.
* means the given filename without suffix.

${txtbld}OPTIONS${txtrst}:
	-f	Data file (Tophat outputted bam or file)
 		${bldred}[NECESSARY]${txtrst}
	-t	File type of input.
 		${bldred}[bam or sam]${txtrst}
	-l	library-type[${bldred}Default fr-secondstrand, acceppt
		fr-firststrand${txtrst}]
	-e	Output reads with wrong mapped flags[${bldred}Default FALSE,
		accept a string to be the name of outputed file${txtrst}]
EOF
}

file=
ltype="fr-secondstrand"
error='FALSE'
itype=''

while getopts "hf:t:l:e:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		t)
			itype=$OPTARG
			;;
		l)
			ltype=$OPTARG
			;;
		e)
			error=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done
if [ -z $file ]; then
	usage
	exit 1
fi

script=${file}".correctTophatOutput.awk"

cat <<EOF >${script}
#!/bin/awk -f

BEGIN{
	OFS="\t";FS="\t";libtype="${ltype}";
	save_discrepancy_to_file="${error}";
	if(save_discrepancy_to_file!="FALSE") system("[ -e " save_discrepancy_to_file " ] && rm " save_discrepancy_to_file);
}
{
    if(\$1 ~ /^@/) print;
    else
    {
        for(i=1;i<=NF;i++) if(\$i!~/^XS/) printf("%s\t",\$i); else XS0=\$i;
        XS1=XS0;
        if(\$2~/^0x/ || \$2~/^[0-9]+\$/){   # FLAG in HEX or Decimal format
            if(libtype=="fr-firststrand") XS1=((and(\$2, 0x10) && and(\$2, 0x40)) || (and(\$2,0x80) && !and(\$2,0x10)))?"XS:A:+":"XS:A:-";
            if(libtype=="fr-secondstrand") XS1=((and(\$2, 0x10) && and(\$2, 0x80)) || (and(\$2,0x40) && !and(\$2,0x10)))?"XS:A:+":"XS:A:-";
        }
        else if(\$2~/^[:alpha:]/){   # FLAG in string
            if(libtype=="fr-firststrand") XS1=(\$2~/r.*1/ || (\$2~/2/ && \$2!~/r/))?"XS:A:+":"XS:A:-";
            if(libtype=="fr-secondstrand") XS1=(\$2~/r.*2/|| (\$2~/1/ && \$2!~/r/))?"XS:A:+":"XS:A:-";
        }
        print XS1;

        if(save_discrepancy_to_file!="FALSE" && XS1!=XS0) print >> save_discrepancy_to_file;
    }
}

EOF

if test $itype == "sam"; then
	output=${file/sam/corrected.sam}
	cat $file | awk -f ${script} >$output 
else
	output_p=${file/bam/corrected.bam}
	samtools view -h ${file} | awk -f ${script} \
	| samtools view -Sb - -o ${output_p}
fi

