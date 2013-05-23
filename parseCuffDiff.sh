#!/bin/bash

#set -x

if test $# -lt 4; then
	echo 1>&2 "Usage $0 parseCuffDiff.py_gene_output \
parseCuffmerge.merged_gtf.py_gname_output \
parseCuffDiff.py_isoform_output \
parseCuffmerge.merged_gtf.py_trname_output"
	exit 1
fi

dir=$(dirname $1)

if test dir == '.'; then
	dir='RNA-Seq/'
else
	dir=${dir}/
fi

gene_exp=$(basename $1)
gene_map=$2

tr_exp=$(basename $3)
tr_map=$4

#---map---change label---------------
#--------------gene--------------------------
awk 'BEGIN{OFS="\t";FS="\t"}ARGIND==1{if(map[$1]=="") map[$1]=$2; \
else map[$1]=map[$1]"&"$2;}\
ARGIND==2{if(FNR==1) print "label",$0; else {\
ne=split(map[$1],name,"&"); for (n in name) print name[n],$0;}}' \
${gene_map} ${dir}${gene_exp} | cut -f 2 --complement >${dir}${gene_exp}.label

awk 'BEGIN{OFS="\t";FS="\t"}{if (FNR==1) \
{output=FILENAME"_g"; print $0 >output; \
output=FILENAME"_x"; print $0 >output; \
output=FILENAME"_u"; print $0 >output; \
output=FILENAME"_c"; print $0 >output; \
output=FILENAME"_o"; print $0 >output; \
output=FILENAME"_s"; print $0 >output; } \
else { a=split($1,b,"__"); cc=b[a]; \
if(cc=="=" || cc=="j") cc="g"; output=FILENAME"_"cc; \
print $0 >>output;}}' ${dir}${gene_exp}.label

for i in `ls ${dir} | grep "${gene_exp}.label_[a-z]$"`; do 
	awk 'BEGIN{OFS="\t";FS="\t"}{if(FNR==1) print $0; else \
	{sum=0;for(i=2;i<=NF;i++) sum+=$i; if (sum>0) print $0;}}' \
	${dir}/$i >${dir}/$i.all_expr.fpkm
done


#--------------isoform----------------------------------
awk 'BEGIN{OFS="\t";FS="\t"}ARGIND==1{if(map[$1]=="") map[$1]=$2; \
else map[$1]=map[$1]"&"$2;}\
ARGIND==2{if(FNR==1) print "label",$0; else {\
ne=split(map[$1],name,"&"); for (n in name) print name[n],$0;}}' \
${tr_map} ${dir}${tr_exp} | cut -f2 --complement >${dir}${tr_exp}.label

awk 'BEGIN{OFS="\t";FS="\t"}{if (FNR==1) \
{output=FILENAME"_a"; print $0 >output; \
output=FILENAME"_j"; print $0 >output; \
output=FILENAME"_x"; print $0 >output; \
output=FILENAME"_u"; print $0 >output; \
output=FILENAME"_c"; print $0 >output; \
output=FILENAME"_o"; print $0 >output; \
output=FILENAME"_s"; print $0 >output; } \
else { a=split($1,b,"__"); cc=b[a]; ; \
if(cc=="=") cc="a"; output=FILENAME"_"cc; \
print $0 >>output;}}' ${dir}${tr_exp}.label

for i in `ls ${dir} | grep "${tr_exp}.label_[a-z]$"`; do 
	awk 'BEGIN{OFS="\t";FS="\t"}{if(FNR==1) print $0; else \
	{sum=0;for(i=2;i<=NF;i++) sum+=$i; if (sum>0) print $0;}}' \
	${dir}/$i >${dir}/$i.all_expr.fpkm
done


#---------------boxplot------------------------------
awk 'BEGIN{OFS="\t";FS="\t"}{if(NR==1) print $0,"Set"; 
  else {print $0,"g"}}' ${dir}${gene_exp}.label_g.all_expr.fpkm \
	>${dir}${gene_exp}.variousType.all_expr.forboxplot
  
for i in `ls ${dir} | grep "${tr_exp}.label_[a-z].all_expr.fpkm$"`; do
	cc=${i/${tr_exp}.label_/}
	cc=${cc/.all_expr.fpkm/}
	awk -v cc=${cc} 'BEGIN{OFS="\t";FS="\t"}{if(NR>1){print $0,cc}}' ${dir}$i \
	>>${dir}${gene_exp}.variousType.all_expr.forboxplot
done

boxplot.onefile.sh -f ${dir}${gene_exp}.variousType.all_expr.forboxplot \
-a Set -s TRUE -v "scale_y_log10()" -S 1 -n FALSE \
-L "'g','a','j','x','s','c','o','u'" \
-p "+theme(legend.position=c(0.08,0.8),legend.title=element_blank())"


