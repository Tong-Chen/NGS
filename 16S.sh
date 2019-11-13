#! /bin/bash
#############################################################
# Title: Pipeline of 16S amplicon by Hiseq2500-PE250
# Author: Yong-Xin Liu
# Wechat: yongxinliu
# E-mail: yxliu@genetics.ac.cn
# Website: http://bailab.genetics.ac.cn/
# Date: 3/16/2017
# Version: 1.3
# Enviroment: Ubuntu 16.04x64, qiime 1.9.1, usearch8
# Description: Script for automatic from clean data and mapping file to otu table, taxonomy, alpha and beta raw result
# Requre file list:
# 1. clean data: clean_data/$library_1/2.fq.gz
# 2. mapping file: doc/$library.mappingfile.txt
#############################################################
# Hiseq2500 PE250 of bacterial 16S, data in clean_data/ and each library mapping file in doc/ 

# 0. Set enviroment and format sample name
wd=/mnt/bai/yongxin/ath/jt.terpene.16S/v1.3 # working directory
# For practice, quickly prepare data and mapping file in your own directory
ln -s /mnt/bai/yongxin/ath/jt.terpene.16S/clean_data/
ln -s /mnt/bai/yongxin/ath/jt.terpene.16S/doc/
# Standard pipeline parameter
rdp=/mnt/bai/public/ref/rdp_gold.fa # rdp gold database, fro remove chimera
gg_align=/mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta # greengene bacterial 16S database
gg_seq=/mnt/bai/public/ref/gg_13_8_otus/rep_set/97_otus.fasta
gg_tax=/mnt/bai/public/ref/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
log=result/readme.log # log file for basic statistics
bc1=0 # forword barcode length, default 0
bc2=6 # forword barcode length, default 0
quality=19 # base quality, accurate > 99%; 29 means 99.9%
bt=6 # barcode type, usually length equal barcode 1 add barcode 2
primer5=AACMGGATTAGATACCCKG # 5` primer used for 16S
primer3=GGAAGGTGGGGATGACGT # 3` primer used for 16S, must reverse compliment
min_len=300 # min length, recommend 300 for bacterial 16S and 220 for ITS
thre_count=3000 # sample min count, filter samples less than thre_count
minuniquesize=2 # min count of unique reads
sim=0.97 # similarity of cluster OTU
p=32 # threads number used: 32
tax_per=0.005 # filter OTU percentage > 0.5% for draw taxonomy and phylogenetic tree, 0.1% about 150 OTU is too much to show

# OTU taxonomy and abundance filter parameter
thre=0.001 # threshold of filter low abundance OTU
taxonomy=p__Cyanobacteria,p__Chloroflexi # filter some phylum
result=result_k1-c # result based on filter OTU table

mkdir $wd
cd $wd

mkdir clean_data
mkdir doc
mkdir temp
mkdir result
# Upload clean data to directory clean_data/, mapping file to each library to doc
# manually modify each library name or using rename regexp batch rename
#rename 's/A201701170/T/g;s/_H55J[3G]BCXY_L[12]//g;s/\.clean//g' clean_data/* 
list=`ls clean_data/*|cut -f 2 -d '/'|cut -f 1 -d '_'|uniq|tr "\n" " "` # library name list



# Batch split each library
for library in $list
	do
	echo -e $(date)"\nStart splitting library $library" >> $log

	# 1. Quality control
	#fastqc -t 9 clean_data/$library\_*.fq.gz

	# 2. Merge clean reads
	join_paired_ends.py -f clean_data/$library\_1.fq.gz -r clean_data/$library\_2.fq.gz -m fastq-join -o temp/$library\_join
	echo 'Merged clean reads:' >> $log
	grep -c -P '^\+$' temp/$library\_join/fastqjoin.join.fastq >> $log

	# 3. Split library
	extract_barcodes.py -f temp/$library\_join/fastqjoin.join.fastq \
		-m doc/$library.mappingfile.txt \
		-o temp/$library\_barcode \
		-c barcode_paired_stitched --bc1_len $bc1 --bc2_len $bc2 -a --rev_comp_bc2
	echo 'Oriented reads:' >> $log
	grep -c -P '^\+$' temp/$library\_barcode/reads.fastq >> $log
	split_libraries_fastq.py -i temp/$library\_barcode/reads.fastq \
		-b temp/$library\_barcode/barcodes.fastq \
		-m doc/$library.mappingfile.txt \
		-o temp/$library\_split/ \
		-q $quality --max_bad_run_length 3 --min_per_read_length_fraction 0.75 --max_barcode_errors 0 --barcode_type $bt
	echo 'Splitted library reads:' >> $log
	grep -c '>' temp/$library\_split/seqs.fna >> $log
	cat temp/$library\_split/split_library_log.txt >> $log # summary library

	# 4. Remove adaptor
	cutadapt -g $primer5 -e 0.15 --discard-untrimmed temp/$library\_split/seqs.fna -o temp/$library\_P5.fa
	echo 'Remove 5` primer:' >> $log
	grep -c '>' temp/$library\_P5.fa >> $log
	cutadapt -a $primer3 -e 0.15 --discard-untrimmed -m $min_len temp/$library\_P5.fa -o temp/$library\_P53.fa
	echo "Remove 3\` primer and length less than $min_len nt:" >> $log
	grep -c '>' temp/$library\_P53.fa >> $log
	#grep -v '>' temp/$library\_P53.fa|awk '{print length($0)}'|sort -n|uniq -c # length distribution
	echo -e "\n\n\n" >> $log

	done

# Merge each library, format to usearch
cat temp/*_P53.fa | sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$/;/g' > temp/seqs_usearch.fa
echo -e $(date)"\nFinished splitting all libraries.\nTotal reads of merge each library :" >> $log
grep -c '>' temp/seqs_usearch.fa >> $log



# 5. Cluster OTU by Usearch
usearch8 -derep_fulllength temp/seqs_usearch.fa \
	-fastaout temp/seqs_unique.fa \
	-minuniquesize $minuniquesize -sizeout >> $log 2>&1
echo 'Unique reads:' >> $log
grep -c '>' temp/seqs_unique.fa >> $log
usearch8 -cluster_otus temp/seqs_unique.fa \
	-otus temp/otus.fa \
	-uparseout temp/otus.up -sizein -sizeout >> $log 2>&1
#echo 'Cluster OTU:' >> $log
#grep -c '>' temp/otus.fa >> $log

# 6. Remove chimeras by rdp_gold database
usearch8 -uchime_ref temp/otus.fa  \
	-nonchimeras temp/otus_rdp.fa \
	-uchimeout temp/otus_rdp.uchime -db $rdp -strand plus >> $log 2>&1
echo 'Remove chimeras by rdp_gold database:' >> $log
grep -c '>' temp/otus_rdp.fa >> $log
align_seqs.py -i temp/otus_rdp.fa -t $gg_align -o temp/aligned/
fasta_subtraction.pl -i temp/otus_rdp.fa -d temp/aligned/otus_rdp_failures.fasta -o temp/otus_rdp_align.fa
echo 'Remove non-bac seq by align_seqs.py:' >> $log
grep '>' -c temp/otus_rdp_align.fa >> $log

# 7. Generate representitive sequences(rep seqs) and OTU table, remove low abundance samples
awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' temp/otus_rdp_align.fa > result/rep_seqs.fa
usearch8 -usearch_global temp/seqs_usearch.fa -db result/rep_seqs.fa -uc temp/otu_table.uc -strand plus -id $sim >> $log 2>&1
python /mnt/bai/public/bin/uc2otutab.py temp/otu_table.uc > temp/otu_table_raw.txt 
biom convert -i temp/otu_table_raw.txt -o temp/otu_table_raw.biom --table-type="OTU table" --to-json
echo 'Summary of otu_table_raw:' >> $log
biom summarize-table -i temp/otu_table_raw.biom >> $log
filter_samples_from_otu_table.py -i temp/otu_table_raw.biom -o result/otu_table.biom -n $thre_count
echo 'Summary of otu_table:' >> $log
biom summarize-table -i result/otu_table.biom >> $log
biom summarize-table -i result/otu_table.biom > result/otu_table.sum
biom convert -i result/otu_table.biom -o result/otu_table.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table.txt

# 8. Taxonomy assignment
assign_taxonomy.py -i result/rep_seqs.fa -r $gg_seq -t $gg_tax -m uclust -o result
sed 's/;/\t/g;s/ //g' result/rep_seqs_tax_assignments.txt > result/rep_seqs_tax.txt # format for R read
mv result/rep_seqs_tax_assignments.log temp/rep_seqs_tax_assignments.log
biom add-metadata -i result/otu_table.biom --observation-metadata-fp result/rep_seqs_tax_assignments.txt -o result/otu_table_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy # add taxonomy to biom
biom convert -i result/otu_table_tax.biom -o result/otu_table_tax.txt --to-tsv --header-key taxonomy
summarize_taxa.py -i result/otu_table_tax.biom -o result/sum_taxa # summary each level percentage
rm result/sum_taxa/*.biom
sed -i '/# Const/d;s/#OTU ID.//g' result/sum_taxa/* # format for R read

# 9. Phylogeny tree
clustalo -i result/rep_seqs.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=$p
filter_alignment.py -i temp/rep_seqs_align.fa -o temp/  # rep_seqs_align_pfiltered.fa, only very short conserved region saved
make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fasta -o result/rep_seqs.tree # generate tree by FastTree

# 10. Alpha diversity
rarefaction=`head -n 7 result/otu_table.sum|tail -n 1|cut -f 3 -d ' '|cut -f 1 -d '.'`
single_rarefaction.py -i result/otu_table.biom -o temp/otu_table_rare.biom -d $rarefaction
alpha_diversity.py -i temp/otu_table_rare.biom -o result/alpha.txt -t result/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree

# 11. Beta diversity
normalize_table.py -i result/otu_table.biom -o temp/otu_table_css.biom -a CSS
biom convert -i temp/otu_table_css.biom -o result/otu_table_css.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' result/otu_table_css.txt
beta_diversity.py -i temp/otu_table_css.biom -o result/beta/ -t result/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
sed -i 's/^\t//g' result/beta/*

# 12. Taxonomy tree - GraPhlAn
filter_otus_from_otu_table.py --min_count_fraction $tax_per -i result/otu_table.biom -o temp/tax_otu_table.biom
filter_fasta.py -f result/rep_seqs.fa -o temp/tax_rep_seqs.fa -b temp/tax_otu_table.biom 
echo "Number of OTU abundance > $tax_per :" >> $log
grep -c '>' temp/tax_rep_seqs.fa >> $log
grep '>' temp/tax_rep_seqs.fa|sed 's/>//g' > temp/tax_rep_seqs.id
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs.id|cut -f 2-3|grep 's__'|sed 's/; /\|/g' > temp/tax_full_anno.txt 
echo "Number of OTU abundance > $tax_per with fully annotation :" >> $log
wc -l temp/tax_full_anno.txt >> $log
echo "Number of OTU abundance > $tax_per with fully annotation unique:" >> $log
sort temp/tax_full_anno.txt|cut -f 1|uniq|wc -l >> $log
format_taxonomy2lefse.pl -i temp/tax_full_anno.txt -o temp/tax_lefse.txt 
## order
export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_order.tree --annotation temp/tax_order.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 4 --min_clade_size 1 --min_font_size 5
graphlan_annotate.py --annot temp/tax_order.annot temp/tax_order.tree temp/tax_order.xml
sed -i 's/ref="A:1">o  /ref="A:1">/g' temp/tax_order.xml
graphlan.py --dpi 300 temp/tax_order.xml result/tax_order.pdf --external_legends
mv result/tax_order_legend.pdf temp/ 
## family
export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_family.tree --annotation temp/tax_family.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 5 --min_clade_size 1 --min_font_size 4
graphlan_annotate.py --annot temp/tax_family.annot temp/tax_family.tree temp/tax_family.xml
sed -i 's/ref="A:1">f  /ref="A:1">/g' temp/tax_family.xml
graphlan.py --dpi 300 temp/tax_family.xml result/tax_family.pdf --external_legends
mv result/tax_family_legend.pdf temp/ 
## genus
export2graphlan.py -i temp/tax_lefse.txt --tree temp/tax_genus.tree --annotation temp/tax_genus.annot --most_abundant 100 --abundance_threshold 0 --least_biomarkers 10 --annotations 6 --min_clade_size 1 --min_font_size 3
graphlan_annotate.py --annot temp/tax_genus.annot temp/tax_genus.tree temp/tax_genus.xml
sed -i 's/ref="A:1">g  /ref="A:1">/g' temp/tax_genus.xml
graphlan.py --dpi 300 temp/tax_genus.xml result/tax_genus.pdf --external_legends
mv result/tax_genus_legend.pdf temp/ 

# 13. Phylogenetic tree - ggtree
clustalo -i temp/tax_rep_seqs.fa -o temp/tax_rep_seqs_clus.fa --seqtype=DNA --full --force --threads=$p
make_phylogeny.py -i temp/tax_rep_seqs_clus.fa -o temp/tax_rep_seqs.tree
sed "s/'//g" temp/tax_rep_seqs.tree > result/tax_rep_seqs.tree # remove '
grep '>' temp/tax_rep_seqs_clus.fa|sed 's/>//g' > temp/tax_rep_seqs_clus.id
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax_assignments.txt temp/tax_rep_seqs_clus.id|sed 's/; /\t/g'|cut -f 1-5 |sed 's/p__//g;s/c__//g;s/o__//g' > result/tax_rep_seqs.tax

# R visualization
cd result
Rscript diversity.R # Draw alpha, beta and Constrain PCoA
Rscript taxonomy.R # Draw barplot+error bar, stack plot
Rscript DEOTU.R # Draw volcano, manhattan, heatmap, venn


# Optional 1. Filter taxonomy and low abundance OTU

# QIIME filter OTU table
mkdir $result
cd result
Rscript filter_otu_table.R
cd ..
filter_otus_from_otu_table.py -i result/otu_table_tax.biom -o temp/k1.biom --otu_ids_to_exclude_fp result/otu_id_k1.txt --negate_ids_to_exclude
echo 'Summary of otu_table_k1, one of sample OTU > 0.1%:' >> $log
biom summarize-table -i temp/k1.biom  >> $log
filter_taxa_from_otu_table.py -i temp/k1.biom -o $result/otu_table.biom -n $taxonomy
echo 'Summary of otu_table_k1 remove:'$taxonomy >> $log
biom summarize-table -i $result/otu_table.biom >> $log
biom summarize-table -i $result/otu_table.biom > $result/otu_table.sum
filter_fasta.py -f result/rep_seqs.fa -o $result/rep_seqs.fa -b $result/otu_table.biom
ln $result/otu_table.biom $result/otu_table_tax.biom
summarize_taxa.py -i $result/otu_table_tax.biom -o $result/sum_taxa
rm $result/sum_taxa/*.biom
sed -i '/# Const/d;s/#OTU ID.//g' $result/sum_taxa/*
biom convert -i $result/otu_table.biom -o $result/otu_table.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU ID.//' $result/otu_table.txt 
cut -f 1 $result/otu_table.txt | tail -n+2 > temp/k1_t.id
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR {a[$1]=$0} NR>FNR {print a[$1]}' result/rep_seqs_tax.txt temp/k1_t.id > $result/rep_seqs_tax.txt

# 9. Phylogeny tree
clustalo -i $result/rep_seqs.fa -o temp/rep_seqs_align.fa --seqtype=DNA --full --force --threads=$p
filter_alignment.py -i temp/rep_seqs_align.fa -o temp/  # rep_seqs_align_pfiltered.fa, only very short conserved region saved
make_phylogeny.py -i temp/rep_seqs_align_pfiltered.fasta -o $result/rep_seqs.tree # generate tree by FastTree

# 10. Alpha diversity
rarefaction=`head -n 7 $result/otu_table.sum|tail -n 1|cut -f 3 -d ' '|cut -f 1 -d '.'`
single_rarefaction.py -i $result/otu_table.biom -o temp/otu_table_rare.biom -d $rarefaction
alpha_diversity.py -i temp/otu_table_rare.biom -o $result/alpha.txt -t $result/rep_seqs.tree -m shannon,chao1,observed_otus,PD_whole_tree

# 11. Beta diversity
normalize_table.py -i $result/otu_table.biom -o temp/otu_table_css.biom -a CSS
biom convert -i temp/otu_table_css.biom -o $result/otu_table_css.txt --table-type="OTU table" --to-tsv
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' $result/otu_table_css.txt
beta_diversity.py -i temp/otu_table_css.biom -o $result/beta/ -t $result/rep_seqs.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
sed -i 's/^\t//g' $result/beta/*

# R visualization
cd $result
Rscript diversity.R # Draw alpha, beta and Constrain PCoA
Rscript taxonomy_bar.R # Draw barplot+error bar, stack plot
Rscript DEOTU.R # Draw volcano, manhattan, heatmap, venn
