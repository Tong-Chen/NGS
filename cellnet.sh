#!/bin/bash

#set -x
set -e
set -u

usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to do CellNet analysis.

Description_file:

1. Normally you should not change the header line especially
'description1' which will be used as an inner variable.

2. 'description1' is used to indicate the replicates of each sample,
same name represents replicate

3. 'file_name' is the actual filename.

4. If 'file_name' column is not supplied, 'sample_id' column will
be treated as filenames without suffix.

#--Below is file content------------------------------------
"sample_id","sample_name","description1","file_name"
"JIAO-D3","JIAO-D3","JIAO-D3","JIAO-D3.cel"
"JIAO-II1-D3","JIAO-II1-D3.cel","JIAO-D3","JIAO-II1-D3.cel"
"JIAO-II2-D3","JIAO-II2-D3.cel","JIAO-D3","JIAO-II2-D3.cel"
"KE-7272-D3","KE-7272-D3.cel","KE-D3","KE-7272-D3.cel"
"MSC-1-D2","MSC-1-D2.cel","MSC-D2","MSC-1-D2.cel"
"MSC-2-D2","MSC-2-D2.cel","MSC-D2","MSC-2-D2.cel"
"PLGA-D3","PLGA-D3.cel","PLGA-D3","PLGA-D3.cel"
"PLGA-II1-D3","PLGA-II1-D3.cel","PLGA-D3","PLGA-II1-D3.cel"
#--Above is file content------------------------------------

${txtbld}OPTIONS${txtrst}:
	-f	Sample description file in format descibed above. ${bldred}[NECESSARY]${txtrst}
	-o	Output prefix for results.${bldred}Default cellnet, [NECESSARY]${txtrst}
	-r	The R command you want to use. Since only R varsion
		2.15.1 to 2.15.3 is suitable for usages, when you have
		multiple R interpretors, please specify the right one by
		giving the absolute path of that command or an alias name to
		represent suitable R version.
		${bldred}[Default Rscript, using current version]${txtrst}
	-l	The library path of needed packages. If not installed, the
		program will try to install all of them when -i is TRUE.
		For my usgae default is "~/home/server-project/cellnet/lib"
		${bldred}[NECESSARY]${txtrst}
	-m	The affy microarray platform used in your project. Currently
		only "hgu133plus2" and "mouse4302" are supported.
		${bldred}[NECESSARY]${txtrst}
	-s	Set the species. Currently only 'human' and 'mouse' are
		supported. ${bldred}[NECESSARY]${txtrst}
	-t	Set the expected target cell type.
		${bldred}[Default 'esc', not NECESSARY]${txtrst}
		Currently avaiblae cell or tissues are:
		mouse: esc, ovary, testis, skin, neuron, glia, hspc,
		macrophage, bcell, tcell, fibroblast, wat, muscleSkel, heart,
		kidney, lung, liver, pancreas, si, colon
	-p	The path for CellNet objects containing the classifiers, GRNs.
		For my usgae default is "~/home/server-project/cellnet/". 
		${bldred}[NECESSARY, needs the last slash]${txtrst}
	-c	The name of the column in the sample data table that indicates
		experimental groups or replicates. Normally the header name of
		the third column should be given here.
		${bldred}[NECESSARY, default 'description1' if header line in
		description file is the same as described above]${txtrst}
	-R	Reload last session if saved.
		${bldred}Default FALSE. One san set -R TRUE when running
		without recomputing.${txtrst}
	-e	Execture the generated R script. [${bldred}Default TRUE${txtrst}]
	-i	Install necessary packages [${bldred}Default FALSE${txtrst}]
EOF
}

file=
R="Rscript"
lib="~/home/server-project/cellnet/lib"
microarrayPlatform=
cellnet_grn="~/home/server-project/cellnet/"
expr_grp='description1'
ist='FALSE'
ext='TRUE'
outprefix='cellnet'
species=''
targetCT='esc'
reload_s='FALSE'

while getopts "hf:o:r:l:m:s:t:p:c:R:e:i:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		o)
			outprefix=$OPTARG
			;;
		r)
			R=$OPTARG
			;;
		l)
			lib=$OPTARG
			;;
		m)
			microarrayPlatform=$OPTARG
			;;
		s)
			species=$OPTARG
			;;
		t)
			targetCT=$OPTARG
			;;
		p)
			cellnet_grn=$OPTARG
			;;
		c)
			exp_grp=$OPTARG
			;;
		R)
			reload_s=$OPTARG
			;;
		e)
			ext=$OPTARG
			;;
		i)
			ist=$OPTARG
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

outprefix=${outprefix}".cellnet"

cat <<EOF >${outprefix}.r

session_file <- paste("${outprefix}", "RData", sep=".")

#R2.15.1-R.2.15.3 is needed if you want to use Affy 1.36.1. I did not
#make it running with these versions.

### Setting up R environment for CellNet

### Set your library path so that it points to the correct platform and annotation libraries:

.libPaths("${lib}")


cn_stdout_all<-function
### save pdfs of standard output
(cnObj,
 ### cnRes object
 cnProc,
 ### CellNet object used to produce cnObj
 tfScores,
 ### result of running cn_nis
 fname_prefix
 ### what to put at the front of the pdffile name
){
  ##<<note function to produce a 'standard' output consisting of:
  ##<< (1) classification heatmap
  ##<< (2) starting and target TCT GRN establishments
  ##<< (3) Aberrant TCT GRN establishments
  ##<< (4) all TCT GRN establishments
  findWidth<-function(panels){
    #  panels<-tmp[['npanels']];
    height<-ceiling(panels/4)*3;
    if(panels<4){
      width<-panels*3;
    }
    else{
      width<-12;
    }
    list(height=height, width=width)
  }
  
  stQuery<-cnObj[['stQuery']];
  cttBest<-cnProc[['grnList']];
  dLevel<-cnObj[['dLevelQuery']];
  
print("Classification heatmap")
  myWidth<-nrow(stQuery)*.4;
  myHeight<-length(cttBest)*.3;
  
  #xtime<-format(Sys.time(), "%Y_%b_%d_%H_%M_%S");
  fname<-paste(fname_prefix, ".pdf",sep='');
  tempTitle<-paste("Classification heatmap ", fname_prefix, sep='');
  pdf(fname, width=8.5, height=11);
  cn_hmClass(cnObj, main=tempTitle);
  
  #set this up so that only the training data from the 'grn' cell type is also shown,
  #along with the query samples
  grnNames<-rownames(cnObj[['normScoresQuery']]);
  for(grnName in grnNames){    
    tSamp<-paste(grnName, "_train", sep='');
    print(cn_barplot_grnSing(cnObj, cnProc, grnName, grnName, bOrder=NULL, norm=T));    
	print(cn_plotnis(tfScores[[grnName]]));
  }
#   print("plot transcriptional regulator scores")
#   print("for now, this only looks at the target cell/tissue type and is a heatmap")
##  tfS<-tfScores[which(tfScores\$grn==targetType),];
##  tfS<-tfS[,1:(ncol(tfS)-1)];
  dev.off();
  fname;
  ### return the file name of the plot pdf file.
}

print("load cellnet package")

library("cellnetr")

print("load appropriate gene annotation sources")

if ("${species}" == "mouse"){
	print("Process mouse data.")
	library("org.Mm.eg.db")
	# mouse4302
	library("${microarrayPlatform}cdf");
	library("${microarrayPlatform}.db");

	# mogene10
	library("mogene10sttranscriptcluster.db");
	library("mogene10stv1cdf");
}else if ("${species}" == "human"){

	print("Process human data.")
	library("org.Hs.eg.db")
	# hgu133plus2
	library("${microarrayPlatform}.db");
	library("${microarrayPlatform}cdf");

	# hugene10
	library("hugene10sttranscriptcluster.db");
	library("hugene10stv1cdf");
}

#set up path the CellNet objects containing the classifiers, GRNs, etc")
#Edit this so that it points to the path where you saved the CellNet objects ")

path_CN_obj <- "${cellnet_grn}"; 

myPlatform <- "${microarrayPlatform}" 

#name of the column in the sample data table that indicates experimental groups 
#or replicates. Set to "", or "sample_name" if there are no replicates

cName <- "${expr_grp}"

#target tissue or cell type.

targetCT<-"${targetCT}"


if (${reload_s} && file.exists(session_file)){
	print("Reload last saved session!")
	load(session_file)
} else {

### Load data and run CellNet

print("Load sample table and fix the data file names ")
#replace 'sampleTab.csv' with your sample table name. See [Step 2 here](http://dev.cellnet.hms.harvard.edu/run/) for a description of the sample table format and an example table.
stQuery<-expr_readSampTab("cellnet.config.csv");
stQuery<-geo_fixNames(stQuery);

print("Load and normalize expression data")
expQuery<-Norm_cleanPropRaw(stQuery, myPlatform);

print("Load right CellNet object")
cnObjName<-switch(myPlatform,
                  mouse4302 = paste(path_CN_obj, "cnProc_mouse4302_062414.R",sep=''),
                  mogene10stv1 = paste(path_CN_obj, "cnProc_mogene_062414.R",sep=''),
                  hgu133plus2 = paste(path_CN_obj, "cnProc_Hg1332_062414.R",sep=''),
                  hugene10stv1 = paste(path_CN_obj, "cnProc_Hugene_062414.R",sep=''));

cnProc<-utils_loadObject(cnObjName);

print("Run CellNet")
cnObj<-cn_apply(expQuery, stQuery, cnProc, dLevelQuery=cName);

} #END reload 

print("Score transcription factors")

grnNames<-rownames(cnObj[['normScoresQuery']]);
targetCT <- grnNames[1]

tfScores<-cn_nis_all(cnObj, cnProc, targetCT);

for(grnName in grnNames) {
	tmp_tfScores<-cn_nis_all(cnObj, cnProc, grnName);
	tfScores[[grnName]] <- tmp_tfScores[[grnName]]
}

print("Visualize results")

cn_stdout_all(cnObj, cnProc, tfScores, "${outprefix}")

print("Output values")

file_expr <- paste("${outprefix}", ".expr.txt", sep="")
write.table(cnObj\$expQuery, file=file_expr, sep="\t", col.names=NA,
row.names=T, quote=F)

file_match_score <- paste("${outprefix}", ".similarity.txt", sep="")
write.table(cnObj\$classRes, file=file_match_score, sep="\t", col.names=NA,
row.names=T, quote=F)

file_grn_score <- paste("${outprefix}", ".GRNscore.txt", sep="")
write.table(cnObj\$normScoreQuery, file=file_grn_score, sep="\t", col.names=NA,
row.names=T, quote=F)

lnames <- names(tfScores)

for(tname in lnames){
	file_netwrok_influence_score <- paste("${outprefix}",tname,"networkInfluScore.txt", sep=".")
	write.table(tfScores[[tname]], file=file_netwrok_influence_score, sep="\t", 
	col.names=NA, row.names=T, quote=F)
}

save.image(file=session_file)

EOF

if [ "${ext}" = 'TRUE' ];then
	${R} ${outprefix}.r
	if [ "$?" == '0' ]; then /bin/rm -f ${outprefix}.r; fi
fi

#########Depleted R codes-------------------
# Gen  expression in the training data set
#pdf(file="mp_rainbowPlot.pdf")
#mp_rainbowPlot(cnProc[['expTrain']],cnProc[['stTrain']],"POU5F1B", dLevel="description1")
#mp_rainbowPlot(cnProc[['expTrain']],cnProc[['stTrain']],"POU5F1")
#dev.off()


#print("Classification heatmap")
#pdf(file="cn_hmClass.pdf")
#cn_hmClass(cnObj);
#dev.off()

#print("Gene regulatory network status of starting cell type (esc) GRN")
#pdf(file="cn_barplot_grnSing.InitialESC.pdf")
#cn_barplot_grnSing(cnObj, cnProc, "esc", c("esc"), bOrder=NULL, norm=T);
#dev.off()

#print("Gene regulatory network status of target cell type (hspc) GRN")
#pdf(file="cn_barplot_grnSing.TargetHSPC.pdf")
#cn_barplot_grnSing(cnObj, cnProc, "fibroblast", c("fibroblast"), bOrder=NULL, norm=T);
#dev.off()

#print("Network influence score of HSPC GRN transcriptional regulators.")
#print("factor heatmap")
#pdf(file="cn_plotnis.pdf")
#cn_plotnis(tfScores, limit=50);
#dev.off()
