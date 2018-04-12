#!/bin/bash

# N cpus executing
ncpu=$1
resolution=$2

i_init=$3
i_end=$4

g_init=$5
g_end=$6

s_init=$7
s_end=$8

echo $*

#N files making up a ca
nfile=1

limit=10000

for (( i=$i_init; i<$i_end; i++)) 
do
i_num=`printf %02d $i`

for (( g=$g_init; g<$g_end; g++)) 
do
g_num=`printf %02d $g`

dir='files/'$resolution'/'$i_num'_'$g_num

if [ -d "$dir" ]; then
#size=`du $dir | cut -f1`
#echo $dir' size: ' ${size}

#if [ "$size" -lt $limit ]; then

#echo rm -rf $dir

#else

echo 'Do a mtree.' $size ./exec_tree_zoom.sh $ncpu $nfile $i $g $s_init $s_end $resolution
./exec_tree_zoom.sh $ncpu $nfile $i $g $s_init $s_end $resolution
#./exec_tree_zoom.sh $ncpu $nfile $i $g 45 50
#./exec_tree_zoom.sh $ncpu $nfile $i $g 40 45
#./exec_tree_zoom.sh $ncpu $nfile $i $g 35 40
#./exec_tree_zoom.sh $ncpu $nfile $i $g 30 35
#./exec_tree_zoom.sh $ncpu $nfile $i $g 25 30
#./exec_tree_zoom.sh $ncpu $nfile $i $g 20 25
#./exec_tree_zoom.sh $ncpu $nfile $i $g 15 20
#./exec_tree_zoom.sh $ncpu $nfile $i $g 10 15
#./exec_tree_zoom.sh $ncpu $nfile $i $g 0 10

#fi	# DO MTREE
fi	# DIRECTORY EXIXST


done
done

