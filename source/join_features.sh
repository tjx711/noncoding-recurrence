#!/bin/bash

##
## This program is used to merge all features data extracted by the scripts:
##   noncoding_extract_features.py //<= predictors
##   noncoding_mut_stats.py        //<= response variable (frequency)
## into a big table which will be loaded into a R data frame. The columns are
## indicated by the array variable $feacols below.
##
## Copyright @ Ramsey Lab, http://lab.saramsey.org
## 03/2016, initial version, Tanjin Xu, tjx711@gmail.com
##


datdir="/home/tanjinxu/Project/noncoding/result"
output="$datdir/cosmic.noncoding.mut.reg331.allsets.tsv"

## dependence/response variable
# frequency/counts within 100bps window by sample freq from COSMIC
samplefreq_f="cosmic_noncoding_mut_reg224.recur.100bp.samplefreq.tsv"
# frequency/counts with 100bps window by number of mutations w/o sample freq
rawfreq_f="cosmic_noncoding_mut_reg330.recur.100bp.rawcounts.tsv"

## observations/features/predictors
features=("recur" \
          "tfbs" \
          "atf3" \
          "cebpb" \
          "cebpd" \
          "creb1" \
          "egr1" \
          "ets1" \
          "maff" \
          "dhs" \
          "gerp" \
          "tss" \
          "gcc" \
          "samplefreq" \
          "rawcounts") 

## headline of the output file
feacols=("id" \
         "category" \
         "tfbs_count" \
         "tfbs_max_score" \
         "tfbs_avg_score" \
         "atf3_count" \
         "atf3_max_score" \
         "atf3_avg_score" \
         "cebpb_count" \
         "cebpb_max_score" \
         "cebpb_avg_score" \
         "cebpd_count" \
         "cebpd_max_score" \
         "cebpd_avg_score" \
         "creb1_count" \
         "creb1_max_score" \
         "creb1_avg_score" \
         "egr1_count" \
         "egr1_max_score" \
         "egr1_avg_score" \
         "ets1_count" \
         "ets1_max_score" \
         "ets1_avg_score" \
         "maff_count" \
         "maff_max_score" \
         "maff_avg_score" \
         "dhs_src_count" \
         "dhs_max_score" \
         "dhs_avg_score" \
         "gerp_score" \
         "tss_dist" \
         "gc_percent" \
         "samplefreq" \
         "rawcounts")

headline=""
for feacol in ${feacols[*]}
do
  if [ $feacol != "rawcounts" ]; then
    headline="$headline""$feacol""\t"
  else
    headline="$headline""$feacol"
  fi
done

echo -e $headline > $output

mv $datdir/result.tsv $datdir/result.tsv.bak

## Catenate all the features together including the responses
for fea in ${features[*]}
do
  #printf "%s\n" $fea
  infile=`find $datdir -name cosmic_noncoding_mut*."$fea"*.tsv`
  infarr=($infile)
  fcnt=${#infarr[*]}

  if [ "$fcnt" -gt 1 ]; then
    echo "# of files found for feature \"$fea\": $fcnt"
    echo $infile
    exit 1
  fi

  if [ -z $infile ]; then
    echo "no file found for feature: \"$fea\""
    exit 1
  fi

  echo $infile
  if [ $fea == "recur" ]; then
    awk -F"\t" '{print $1"-"$2"\t"$4}' $infile | sort -k1,1 > "$datdir/result.tsv"
  elif [ $fea == "dhs" ]; then
    awk -F"\t" '{print $1"-"$2"\t"$4"\t"$5"\t"$6}' $infile|sort -k1,1 > "$datdir/tmp2"
  elif [ $fea == "tss" ] || [ $fea == "gcc" ]; then
    awk -F"\t" '{print $1"-"$2"\t"$3}' $infile|sort -k1,1 > "$datdir/tmp2"
  elif [ $fea == "gerp" ]; then
    awk -F"\t" '{print $1"\t"$4}' $infile|sort -k1,1 > "$datdir/tmp2"
  elif [ $fea == "samplefreq" ] || [ $fea == "rawcounts" ]; then
    awk -F"\t" '{print $1"-"$2"\t"$3}' $infile|sort -k1,1 > "$datdir/tmp2"
  else
    awk -F"\t" '{print $1"-"$2"\t"$3"\t"$5"\t"$6}' $infile|sort -k1,1 > "$datdir/tmp2"
  fi

  if [ $fea != "recur" ]; then
    join -t$'\t' -a1 -1 1 -2 1 -o auto -e "0" "$datdir/result.tsv" "$datdir/tmp2" > "$datdir/result2.tsv"
    mv "$datdir/result2.tsv" "$datdir/result.tsv"
    #rm -rf "$datdir/tmp2"
  fi

  fsize=`wc -l $datdir/result.tsv`
  ncols=`awk -F"\t" '{print NF}' $datdir/result.tsv|head -n 1`
  echo "feature:$fea,$ncols,$fsize"
done

# sort the merged data by chromosome
for chr in {1..22} {X,Y}
do 
  grep "^$chr-" $datdir/result.tsv >> $output
done

