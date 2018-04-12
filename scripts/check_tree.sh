#!/bin/bash

resolution=$1

i_init=$2
i_end=$3

g_init=$4
g_end=$5

s_init=$6
s_end=$7

for (( i=$i_init; i<$i_end; i++)) 
do
i_num=`printf %02d $i`

for (( g=$g_init; g<$g_end; g++)) 
do
g_num=`printf %02d $g`

dir='output/'$resolution'/'$i_num'_'$g_num

if [ -d "$dir" ]; then
echo $dir

for (( s=$s_init; s<$s_end; s++ ))
do
s_num=`printf %03d $s`
s_nump=`printf %03d $(( $s + 1 ))`

file=`ls $dir/bothways_${s_nump}_${s_num}.0000_mtree`

done

fi

done
done

