#!/usr/bin/env bash

set -ue

if [ $# -lt 5 ]
then
  echo "
  SEACR: Sparse Enrichment Analysis for CUT&RUN
  
  Usage: bash SEACR_1.3.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDR threshold>] ["norm" | "non"] ["relaxed" | "stringent"] output prefix
  
  Description of input fields:
  
  Field 1: Target data bedgraph file in UCSC bedgraph format (https://genome.ucsc.edu/goldenpath/help/bedgraph.html) that omits regions containing 0 signal.
  
  Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling. Alternatively, a numeric threshold n between 0 and 1 returns the top n fraction of peaks based on total signal within peaks.
  
  Field 3: “norm” denotes normalization of control to target data, “non” skips this behavior. "norm" is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).
    
  Field 4: “relaxed” uses a total signal threshold between the knee and peak of the total signal curve, and corresponds to the “relaxed” mode described in the text, whereas “stringent” uses the peak of the curve, and corresponds to “stringent” mode.
  
  Field 5: Output prefix
  
  Output file:
  <output prefix>.auc.threshold.merge.bed (Bed file of enriched regions)
  
  Output data structure: 
  
  <chr>  <start>  <end>  <AUC>  <max signal>  <max signal region>
  
  Description of output fields:
  Field 1: Chromosome
  
  Field 2: Start coordinate
  
  Field 3: End coordinate
  
  Field 4: Total signal contained within denoted coordinates
  
  Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates
  
  Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal
  
  Examples:
  bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output
  Calls enriched regions in target data using normalized IgG control track with stringent threshold
  
  bash SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output
  Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold
  bash SEACR_1.3.sh target.bedgraph 0.01 non stringent output
  Calls enriched regions in target data by selecting the top 1% of regions by area under the curve (AUC)
  "
  exit 1
fi

filter_auc_bed() {
  # fdr, threshold s, threshold g, password
  echo "Empirical false discovery rate = $1"
  awk -v value=$2 -v value2=$3 '$4 > value && $7 > value2 {print $0}' $4.auc.bed | \
     cut -f 1,2,3,4,5,6 > $4.auc.threshold.bed
}

password=`head /dev/urandom | LC_CTYPE=C tr -dc A-Za-z0-9 | head -c 13; echo ''`
password2=`head /dev/urandom | LC_CTYPE=C tr -dc A-Za-z0-9 | head -c 13; echo ''`

exp=`basename $1`

if [[ -f $2 ]]
then
  echo "Calling enriched regions with control file"
  ctrl=`basename $2`
elif [[ $2 =~ ^(1(\.0*)?|0?\.[0-9]+)$ ]]
then
  echo "Calling enriched regions without control file"
else
  echo "$2 is not a file or a number between 0 and 1"
  exit 1
fi

norm=`echo $3`

if [[ $norm == "norm" ]]
then
  echo "Normalizing control to experimental bedgraph"
elif [[ $norm == "non" ]]
  then
  echo "Proceeding without normalization of control to experimental bedgraph"
else
  echo "Must specify \"norm\" for normalized or \"non\" for non-normalized data processing in third input"
  exit 1
fi

height=`echo $4`

if [[ $height == "relaxed" ]]
then
  echo "Using relaxed threshold"
elif [[ $height == "stringent" ]]
  then
  echo "Using stringent threshold"
else
  echo "Must specify \"stringent\" or \"relaxed\" in fourth input"
  exit 1
fi

echo "Creating experimental AUC file: $(date)"
python3 bdg_to_auc.py $1 > $password.auc.bed

if [[ -f $2 ]]
then
  echo "Creating control AUC file: $(date)"
  python3 bdg_to_auc.py $2 > $password2.auc.bed
fi

echo "Calculating optimal AUC threshold: $(date)"

if [[ -f $2 ]] && [[ $norm == "norm" ]]
then
  echo "Calculating threshold using normalized control: $(date)"
  Rscript SEACR_1.3.new.R --exp=$password.auc.bed --ctrl=$password2.auc.bed --norm=yes --output=$password
elif [[ -f $2 ]]
then
  echo "Calculating threshold using non-normalized control: $(date)"
  Rscript SEACR_1.3.new.R --exp=$password.auc.bed --ctrl=$password2.auc.bed --norm=no --output=$password
else
  echo "Using user-provided threshold: $(date)"
  Rscript SEACR_1.3.new.R --exp=$password.auc.bed --ctrl=$2 --norm=no --output=$password
fi

source ${password}_var.sh
echo "Creating thresholded feature file: $(date)"

if [[ $height == "relaxed" ]]
then
  filter_auc_bed $SEACR_FDR_R $SEACR_THRESH_R $SEACR_THRESH_G $password
else
  filter_auc_bed $SEACR_FDR_S $SEACR_THRESH_S $SEACR_THRESH_G $password
fi

if [[ -f $2 ]]
then
  if [[ $norm == "norm" ]] #If normalizing, multiply control bedgraph by normalization constant
  then
    awk -v mult=$SEACR_NORM 'BEGIN{OFS="\t"}; {$4=$4*mult; print $0}' $password2.auc.bed | cut -f 1,2,3,4,5,6 > $password2.auc2.bed
    mv $password2.auc2.bed $password2.auc.bed
  fi
  awk -v value=$SEACR_THRESH_S '$4 > value {print $0}' $password2.auc.bed > $password2.auc.threshold.bed
fi

echo "Merging nearby features and eliminating control-enriched features: $(date)"

# module load bedtools ## For use on cluster
mean=`awk '{s+=$3-$2; t++}END{print s/(t*10)}' $password.auc.threshold.bed`

if [[ -f $2 ]]
then
  awk -v value=$mean 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6; s++}else{if(chr==$1 && $2 < stop+value){stop=$3; auc=auc+$4; if($5 > max){max=$5; coord=$6}else if($5==max){split(coord,t,"-"); split($6,u,"-"); coord=t[1]"-"u[2]}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6}}}' $password.auc.threshold.bed | bedtools intersect -wa -v -a - -b $password2.auc.threshold.bed > $5.auc.threshold.merge.bed  
else
  awk -v value=$mean 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6; s++}else{if(chr==$1 && $2 < stop+value){stop=$3; auc=auc+$4; if($5 > max){max=$5; coord=$6}else if($5==max){split(coord,t,"-"); split($6,u,"-"); coord=t[1]"-"u[2]}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6}}}' $password.auc.threshold.bed > $5.auc.threshold.merge.bed
fi

if [[ $height == "relaxed" ]]
then
  cat $5.auc.threshold.merge.bed > $5.relaxed.bed
else
  cat $5.auc.threshold.merge.bed > $5.stringent.bed
fi

echo "Removing temporary files: $(date)"
rm ${password}_var.sh

#rm $password.auc.bed
#rm $password.auc
#rm $password.auc.threshold.bed
rm $5.auc.threshold.merge.bed
if [[ -f $2 ]]
then
  rm $password2.auc.bed
  rm $password2.auc
  rm $password2.auc.threshold.bed
fi

echo "Done: $(date)"
