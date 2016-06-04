#!/bin/bash

set -x;
datdir=/home/tanjinxu/Project/noncoding/data
rfl1=$datdir/cosmic_noncoding_mut.tsv
rfl2=$datdir/cosmic_noncoding_mut_reg201.recur.tsv

tab="$(printf \t)"
rm -rf mut.tmp.out
for chr in {1..22} {X,Y}
do
  echo $chr
  awk -F"\t" -v var1="$chr" '{if ($1 == var1) print}' $rfl1 | sort -k2,2 > tmp1
  awk -F"\t" -v var2="$chr" '{if ($1 == var2) print}' $rfl2 | sort -k2,2 > tmp2
  join -1 2 -2 2 -o 1.1,1.2,2.3,2.4,2.5,1.3,1.4 tmp2 tmp1 |sort -n -k2,2 >>mut.tmp.out
done
