#!/bin/bash

resolution=$1

cpus=$2

i_init=$3
i_end=$4

g_init=$5
g_end=$6

s_init=$7
s_end=$8

for (( i=$i_init; i<$i_end; i++)) 
do
i_num=`printf %02d $i`

for (( g=$g_init; g<$g_end; g++)) 
do
g_num=`printf %02d $g`

dir='/home/uam37/uam37526/CLUES/merger_tree/output/'$resolution'/'$i_num'_'$g_num

if [ -d "$dir" ]; then
echo $dir

for (( s=$s_init; s<$s_end; s++ ))
do
s_num=`printf %03d $s`
s_nump=`printf %03d $(( $s + 1 ))`

file=`ls $dir/bothways_${s_nump}_${s_num}.0000_mtree`
#echo $file
if [ -e "$file" ]; then
#do nothing
mmm=0
else

echo 'Fix the trees'
echo ./loop_mtree_sergey.sh $cpus $resolution $i $(( i+1 )) $g $(( g+1 )) $s $(( s+1 ))
#timeout 18s 
timeout 20s ./loop_mtree_sergey.sh $cpus $resolution $i $(( i+1 )) $g $(( g+1 )) $s $(( s+1 ))
#./loop_mtree_sergey.sh $cpus $resolution $i $(( i+1 )) $g $(( g+1 )) $s $(( s+1 ))

fi

done

fi

done
done

