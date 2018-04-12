#!/bin/bash

# N cpus executing
ncpu=$1

# Lowest level of refinement
resolution=$2

i_init=$3
i_end=$4

g_init=$5
g_end=$6

s_init=$7
s_end=$8

echo 'Input params: '$*

#N files making up a catalog
nfile=1

limit=1000

for (( i=$i_init; i<$i_end; i++)) 
do
i_num=`printf %02d $i`

for (( g=$g_init; g<$g_end; g++)) 
do
g_num=`printf %02d $g`

dir='files/'$resolution'/'$i_num'_'$g_num

echo $dir 

if [ -d "$dir" ]; then
size=`du $dir | cut -f1`
echo $dir' size: ' ${size}

if [ "$size" -lt $limit ]; then

echo rm -rf $dir

else

echo 'Do a mtree.' $size
echo ./exec_tree_zoom.sh $ncpu $nfile $i $g $s_init $s_end $resolution
./exec_tree_zoom.sh $ncpu $nfile $i $g $s_init $s_end $resolution 

fi	# DO MTREE
else

echo 'Directory ' $dir ' does not exist.'
fi	# DIRECTORY EXIXST

done
done
