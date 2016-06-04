#!/bin/bash

wdir="/home/tanjinxu/Project/noncoding/plots"
datpath="$wdir/data"
tmpdir="$wdir/tmp"

metrics=("test.mse" \
         "test.mean_ae|test.muae" \
         "test.median_ae|test.meae" \
         "test.r2" \
         "test.expvar|test.eas")

for metric in ${metrics[*]}
do
 datfiles="$datpath/*tsv"

 metstr=""
 if [ $metric == "test.mean_ae|test.muae" ]; then
  metstr="test.muae"
 elif [ $metric == "test.median_ae|test.meae" ]; then
  metstr="test.meae"
 elif [ $metric == "test.expvar|test.eas" ]; then
  metstr="test.expvar"
 else
  metstr=$metric
 fi

 resfile="$wdir/result/nc.res.result.$metstr"
 egrep $metric $datfiles|egrep -v "poisson|nb|quasi"|awk -F":" '{print $1"\t"$2}'|cut -f3-8 > "$tmpdir/tmp.col2"
 egrep $metric $datfiles|egrep -v "poisson|nb|quasi"|awk -F":" '{print $1"\t"$2}'|cut -f1 > "$tmpdir/tmp.col1"
 awk -F"." '{print $5"."$6"\t"$8}' "$tmpdir/tmp.col1" > "$tmpdir/tmp.col3"
 paste "$tmpdir/tmp.col3"  "$tmpdir/tmp.col2" |sort -k1,1 -k2,2 > "$tmpdir/result.sklearn"
 awk -F"\t" '{avg = ($4+$5+$6+$7+$8)/5; print $1"\t"$2"\t"$3"\t"avg}' $tmpdir/result.sklearn > $tmpdir/result
 

 egrep $metric $datfiles|grep -v log|egrep "poisson|nb|quasi"|awk -F":" '{print $1"\t"$2}'|grep -v shuffled|cut -f3 >"$tmpdir/tmp2.pos.nb"
 egrep $metric $datfiles|grep -v log|egrep "poisson|nb|quasi"|awk -F":" '{print $1"\t"$2}'|grep shuffled|cut -f3 > "$tmpdir/tmp2.pos.nb.shuff"
 egrep $metric $datfiles|grep -v log|egrep "poisson|nb|quasi"|awk -F":" '{print $1"\t"$2}'|grep shuffled|cut -f1|awk -F"." '{print $5"."$6"\t"$8}' >"$tmpdir/tmp1.pos.nb"

 paste "$tmpdir/tmp1.pos.nb" "$tmpdir/tmp2.pos.nb" > $tmpdir/result.nb.pos
 paste $tmpdir/result.nb.pos $tmpdir/tmp2.pos.nb.shuff > $tmpdir/tmp
 mv $tmpdir/tmp $tmpdir/result.nb.pos

 cat $tmpdir/result.nb.pos >> $tmpdir/result
 sort -k1,1 -k2,2 $tmpdir/result > $resfile
done         
