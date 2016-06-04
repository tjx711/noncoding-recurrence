#!/bin/bash

wdir=/home/tanjinxu/winsys/research/tanjin/mutations/noncoding
fanno=$wdir/ensembl_hg19_anno117.tsv
chranno_pre=$wdir/ensembl_hg19_anno

set -x;

for i in {1..22} {X,Y}
do
  #echo $i
  chranno_fil=$chranno_pre"_chr"$i".tsv"
  echo $chranno_fil
  if [ -f $chranno_fil ]; then
    echo "the file already exits: $chranno_fil"
    rm -rf $chranno_fil
  fi

  touch $chranno_fil
  #head -n 1 $fanno > $chranno_fil
  #awk -v chr=$i '{if ($3 == chr) print}' $fanno | sort -n -k4,4 -k6,6 >> $chranno_fil
  cut -f1,3,4,5,6,7-11 $fanno |head -n 1 > $chranno_fil
  grep -v -i patch $fanno|cut -f1,3,4,5,6,7-11 |awk -v chr=$i -F"\t" '{if ($2 == chr) print}'|sort -n -k3,3|uniq >> $chranno_fil
done
