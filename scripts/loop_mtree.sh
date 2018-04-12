#!/bin/bash

i_init=0
i_end=80

g_init=0
g_end=10

s_init=0
s_end=5

# N cpus executing
ncpu=$1

#N files making up a ca
nfile=1

limit=1000

resolution='1024lgf'

echo $*

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

#echo rm -rf $dir

echo 'StoCazzo'

else

echo 'Do a mtree.' $size
echo ./exec_tree_zoom.sh $ncpu $nfile $i $g $s_init $s_end $resolution
timeout 5s ./exec_tree_zoom.sh $ncpu $nfile $i $g $s_init $s_end $resolution 

fi	# DO MTREE
fi	# DIRECTORY EXIXST


done
done

