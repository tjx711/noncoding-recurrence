#!/bin/sh

mutfile=data/cosmic_noncoding_mut_reg201.recur.tsv
tfsfile=data/wgEncodeRegTfbsClusteredV3.bed
#dhsfile=wgEncodeRegDnaseClusteredV3.bed
#outfile=cosmic_noncoding_mut_reg205.features.clustered.tsv
dhsfile=data/wgEncodeAwgDnaseMasterSites.bed
#outfile=cosmic_noncoding_mut_reg205.features.master.tsv
#dhsfile=wgEncodeAwgDnaseUniformSites.bed

#etsfile=data/wgEncodeAwgTfbsAtf3.bed
#etsfile=data/wgEncodeAwgTfbsCebpb.bed
#etsfile=data/wgEncodeAwgTfbsCebpd.bed
#etsfile=data/wgEncodeAwgTfbsCreb1.bed
etsfile=data/wgEncodeAwgTfbsEgr1.bed
#etsfile=data/wgEncodeAwgTfbsEts1.bed
#etsfile=data/wgEncodeAwgTfbsMaff.bed

#tssfile=data/hg19_transcripts_tss_d0218.tsv
gccfile=/home/tanjinxu/Project/noncoding/mutations/noncoding/features/gcp/hg19.gc5Base.221.txt
outfile=result/cosmic_noncoding_mut_reg222.gcc.tsv

python noncoding_extract_features.py \
       -gcc \
       $mutfile \
       $gccfile \
       $outfile \
