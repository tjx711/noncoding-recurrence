#!/bin/bash

gerp_data=/home/tanjinxu/Project/noncoding/mutations/noncoding/features/gerp/All_hg19_RS.bw
mut_data=/home/tanjinxu/Project/noncoding/data/cosmic_noncoding_reg_mut.bed
res=~/Project/noncoding/result/cosmic_noncoding_mut_reg218.gerp.tsv

cd /home/tanjinxu/Project/noncoding

nohup ./bigWigAverageOverBed $gerp_data $mut_data $res >gerp_label.out 2>&1 &
