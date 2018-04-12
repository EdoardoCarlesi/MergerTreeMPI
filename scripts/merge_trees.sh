#!/bin/bash

base_tree='/home/uam37/uam37526/scratch/CLUES/merger_tree/output/'
tree_start=0
tree_end=2

ginn=0
run_start=0
run_end=1

root_tree='bothways'
#root_tree='mtree'

ginn_dir=`printf %02d ${ginn}`

for (( irun=$run_start; irun<$run_end; irun++ ))
do
run_dir=`printf %02d ${irun}`

tree_dir=$base_tree'/512_'$ginn_dir'/'$run_dir'/'

for (( itree=$tree_start; itree<$tree_end; itree++ ))
do 

iplustree=`expr ${itree} + 1`

ip=`printf %03d ${iplustree}`
im=`printf %03d ${itree}`

out_mtree=$tree_dir'/mtmerge_'$ip'_'$im'_mtree'
out_mtree_idx=$tree_dir'/mtmerge_'$ip'_'$im'_mtree_idx'

echo cat ${tree_dir}$root_tree*${ip}'_'${im}*'mtree' > $out_mtree
cat ${tree_dir}$root_tree*${ip}'_'${im}*'mtree' > $out_mtree

echo cat ${tree_dir}$root_tree*${ip}'_'${im}*'mtree_idx' > $out_mtree_idx
cat ${tree_dir}$root_tree*${ip}'_'${im}*'mtree_idx' > $out_mtree_idx

done
done
