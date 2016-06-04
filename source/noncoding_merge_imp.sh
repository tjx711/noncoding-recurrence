#!/bin/bash

wdir="/home/tanjinxu/Project/noncoding/plots"
datpath="$wdir/data"
tmpdir="$wdir/tmp"


impfiles=`find $datpath -name "*.imp.tsv" |egrep "nb|poisson"`
impfarr=($impfiles)

for fimp in ${impfarr[*]}
do
  model=`echo $fimp | cut -d"." -f8`
  rate=`echo $fimp | cut -d"." -f5,6`

  scores=`sed -n '1!p' $fimp |cut -f2`
  scorelist=($scores)
  
  #record=$(printf  "%s\t%s\t%s" $model $rate "${scorelist[@]}")
  record=$(printf ",%s" "${scorelist[@]}")
  record=${record:1}
  printf "%s,%s\n" $model $rate >>tmp
  echo $record >> res.imp
  
done
