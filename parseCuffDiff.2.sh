#!/bin/bash

set -x
set -e
set -u

if test $# -lt 8; then
	echo 1>&2 "This has different usage compared with
parseCuffDiff.sh. These two programs are independent."
	echo 1>&2 "Usage $0 gene_exp.diff \
parseCuffmerge.merged_gtf.py_gname_output \
isoform_exp.diff \
parseCuffmerge.merged_gtf.py_trname_output \
gene_prefix isoform_prefix \
least_rpkm least_len"
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
gene_p=$5

tr_exp=$(basename $3)
tr_map=$4
tr_p=$6

least_rpkm=$7
least_len=$8
##---map---change label---------------
##--------------gene--------------------------
echo ">>> Substitute labeled name for ${gene_exp}"
echo ">>>"
awk 'BEGIN{OFS="\t";FS="\t"}ARGIND==1{if(map[$1]=="") map[$1]=$2; \
else map[$1]=map[$1]"&"$2;}\
ARGIND==2{if(FNR==1) print $0; else \
{ne=split(map[$1],name,"&"); for (n in name) {$1=name[n]; print $0;}}}' \
${gene_map} ${dir}${gene_exp} >${dir}${gene_exp}.label

echo ">>> Sort genes with specific labels to separate files."
echo ">>>"
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

echo ">>> Transfer genes into normal table and select genes meets
given conditions. Files with only one line will be skipped."
echo ">>>"
for i in `ls ${dir} | grep "${gene_exp}.label_[a-z]$"`; do 
	if test $(cat ${dir}${i} | wc -l) -gt 1; then
		cc=${i/${gene_exp}.label_/}
		cuffdiff.transfer_diff.py ${dir}$i \
		${gene_p}_${cc} ${least_rpkm} ${least_len}
		cut -f 2-4 --complement \
		${gene_p}_${cc}.expr_${least_rpkm}.len_${least_len} > \
			${gene_p}_${cc}.expr_${least_rpkm}.len_${least_len}.tmp
		/bin/mv -f ${gene_p}_${cc}.expr_${least_rpkm}.len_${least_len}.tmp \
			${gene_p}_${cc}.expr_${least_rpkm}.len_${least_len}
	fi
done

echo ">>> Get expression of genes."
echo ">>>"
gene_gp=${gene_p}_g
gene_p_dir=$(dirname ${gene_gp})
gene_p_base=$(basename ${gene_gp})
if test gene_p_dir == '.'; then
	gene_p_dir='/'
else
	gene_p_dir=${gene_p_dir}/
fi

awk 'BEGIN{OFS="\t";FS="\t"}{if(FNR==1) {for(i=1;i<=NF;i++) {a[i]=$i;if($i~/___/) {out=ARGV[1]"."a[i]; print $0 >out;}}} else {for(i=2;i<=NF;i++) {if ($i=="yes") {out=ARGV[1]"."a[i]; print $0 >>out}}}}' ${gene_gp}.expr_${least_rpkm}.len_${least_len}

echo "Separate genes into files with sample-sample comparison data."
 cuffdiff.systemtial_DE.py $(ls ${gene_p_dir} | grep '___' | grep -P "^${gene_p_base}.expr_${least_rpkm}.len_${least_len}.[^.]*$" | sed "s#^#${gene_p_dir}#")

#--------------isoform----------------------------------

echo ">>> Substitute labeled name for ${tr_exp}"
echo ">>>"
awk 'BEGIN{OFS="\t";FS="\t"}ARGIND==1{if(map[$1]=="") map[$1]=$2; \
else {print "This should never happen."; exit;}}\
ARGIND==2{if(FNR==1) print $0; else \
{$1=map[$1]; print $0;}}' \
${tr_map} ${dir}${tr_exp} >${dir}${tr_exp}.label

echo ">>> Sort out genes into separate files by different labels"
echo ">>>"
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

echo ">>> Transfer isoforms into normal table and select isoforms meets
given conditions. Files with only one line will be skipped."
echo ">>>"
for i in `ls ${dir} | grep "${tr_exp}.label_[a-z]$"`; do 
	if test $(cat ${dir}${i} | wc -l) -gt 1; then
		cc=${i/${tr_exp}.label_/}
		cuffdiff.transfer_diff.py ${dir}$i \
		${tr_p}_${cc} ${least_rpkm} ${least_len}
		cut -f 2-4 --complement \
		${tr_p}_${cc}.expr_${least_rpkm}.len_${least_len} > \
			${tr_p}_${cc}.expr_${least_rpkm}.len_${least_len}.tmp
		/bin/mv -f ${tr_p}_${cc}.expr_${least_rpkm}.len_${least_len}.tmp \
			${tr_p}_${cc}.expr_${least_rpkm}.len_${least_len}
	fi
done

echo ">>> Sort all expressed isoforms"
echo ">>>"
tr_ap=${tr_p}_a

awk 'BEGIN{OFS="\t";FS="\t"}{if(FNR==1) {for(i=1;i<=NF;i++) {a[i]=$i;if($i~/___/) {out=ARGV[1]"."a[i]; print $0 >out;}}} else {for(i=2;i<=NF;i++) {if ($i=="yes") {out=ARGV[1]"."a[i]; print $0 >>out}}}}' ${tr_ap}.expr_${least_rpkm}.len_${least_len}

tr_ap_dir=$(dirname ${tr_ap})

if test tr_ap_dir == '.'; then
	tr_ap_dir='/'
else
	tr_ap_dir=${tr_ap_dir}/
fi
tr_ap_base=$(basename ${tr_ap})

 cuffdiff.systemtial_DE.py $(ls ${tr_ap_dir} | grep '___' | grep -P "^${tr_ap_base}.expr_${least_rpkm}.len_${least_len}.[^.]*$" | sed "s#^#${tr_ap_dir}#")


echo ">>> Sort new isoforms if exists"
echo ">>>"
tr_jp=${tr_p}_j
if test -s ${tr_jp}; then
	#
	#
	awk 'BEGIN{OFS="\t";FS="\t"}{if(FNR==1) {for(i=1;i<=NF;i++) {a[i]=$i;if($i~/___/) {out=ARGV[1]"."a[i]; print $0 >out;}}} else {for(i=2;i<=NF;i++) {if ($i=="yes") {out=ARGV[1]"."a[i]; print $0 >>out}}}}' ${tr_jp}.expr_${least_rpkm}.len_${least_len}

	tr_jp_dir=$(dirname ${tr_jp})

	if test tr_jp_dir == '.'; then
		tr_jp_dir='/'
	else
		tr_jp_dir=${tr_jp_dir}/
	fi
	tr_jp_base=$(basename ${tr_jp})

	 cuffdiff.systemtial_DE.py $(ls ${tr_jp_dir} | grep '___' | grep -P "^${tr_jp_base}.expr_${least_rpkm}.len_${least_len}.[^.]*$" | sed "s#^#${tr_jp_dir}#")
fi
#---------------boxplot----only for three samples--------------------------
#awk 'BEGIN{OFS="\t";FS="\t"}{if(NR==1) print $1,$2,$3,$4,"Set"; \
#	else {print $1,$2,$3,$4,"g"}}' ${gene_gp}.expr_${least_rpkm}.len_${least_len} \
#	>${gene_gp}.expr_${least_rpkm}.len_${least_len}.variousType.all_expr.forboxplot
#  
#
#for i in \
#	`ls ${dir}${tr_ap_dir} | grep \
#	"^${tr_ap_base}_[a-z].expr_${least_rpkm}.len_${least_len}$"`; do
#	cc=${i/${tr_ap_base}_/}
#	cc=${cc/.expr_${least_rpkm}.len_${least_len}/}
#	awk -v cc=${cc} 'BEGIN{OFS="\t";FS="\t"}\
#	{if(NR>1){print $1,$2,$3,$4,cc}}' ${dir}${tr_ap_dir}$i \
#	>>${gene_gp}.expr_${least_rpkm}.len_${least_len}.variousType.all_expr.forboxplot
#done
#
#boxplot.onefile.sh -f \
#${gene_gp}.expr_${least_rpkm}.len_${least_len}.variousType.all_expr.forboxplot \
#-a Set -s TRUE -v "scale_y_log10()" -S 1 -n FALSE \
#-L "'g','a','j','x','s','c','o','u'" \
#-p "+theme(legend.position=c(0.08,0.8),legend.title=element_blank())"


