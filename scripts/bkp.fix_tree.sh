#!/bin/bash

# Merger Tree directories and catalogue path
home_dir='/home/uam37/uam37526/scratch/CLUES/'
#home_dir='/home/edoardo/CLUES/'
mtree=$home_dir'/merger_tree/'
#mtree=$home_dir'/MergerTree/'
mpi_exec='mpiexec'

# N_proc sets the number of mpi tasks
N_proc=5

# Chunks of file
N_read=16

ice_seed=`printf %02d $1`
gin_seed=`printf %02d $2`

seeds=$ice_seed'_'$gin_seed

snap_min=$3
snap_max=$4

echo $N_proc $N_read $seeds $snap_min $snap_max

ls 'output/'$ice_seed'_'$gin_seed'/'*'0000'*'mtree' | wc -l
ls 'output/'$ice_seed'_'$gin_seed'/'*'0000'*'mtree'


# AHF settings
executable=$mtree'/MergerTree'

# Store temporary files
temp_dir=$mtree'/temp/'$ice_seed$gin_seed$snap_min$snap_max

mkdir $temp_dir

# Where to dump the lists of the _halo and _particles files
temp_part=$temp_dir'/list_part_files_'$ice_seed$gin_seed$snap_min$snap_max'.ahf'
temp_halo=$temp_dir'/list_halo_files_'$ice_seed$gin_seed$snap_min$snap_max'.ahf'
temp_out=$temp_dir'/list_files_numbers_'$ice_seed$gin_seed$snap_min$snap_max'.ahf'

# Box size in kpc (needed for the haloes at the box boundaries, for periodic boxes)
BOX=100000
mkdir $mtree'/output/'$seeds'/'

# Output files directory + prefix
mkdir $mtree'/output/'$seeds'/'

output_dir=$mtree'/output/'$seeds'/'
base_catalog=$mtree'/files/'$seeds'/'

output_file=$output_dir'bothways'

edit_temp=0	

rm -rf $temp_part $temp_halo $temp_out
touch $temp_part $temp_halo $temp_out

if [ "$edit_temp" -eq "1" ] ; then
if [ "$2" -eq "1" ] ; then 

ls -r $base_catalog*particles > $temp_part
ls -r $base_catalog*halos > $temp_halo
ls -r $base_catalog*particles | grep -o '_[0-9][0-9][0-9]' > $temp_out

else 

ls -r $base_catalog*0000*particles > $temp_part
ls -r $base_catalog*0000*halos > $temp_halo
ls -r $base_catalog*0000*particles | grep -o '_[0-9][0-9][0-9]' > $temp_out

fi

else 	

for (( isnap=$snap_max; isnap>=$snap_min; isnap-- )) 
do
snap_num=`printf %03d $isnap` 
#echo 'snapshots= ' $snap_num
if [ "$isnap" -eq "$snap_max" ]; then 
ls -r $base_catalog*$snap_num*0000*particles > $temp_part
ls -r $base_catalog*$snap_num*0000*halos > $temp_halo
ls -r $base_catalog*$snap_num*0000*particles | grep -o '_[0-9][0-9][0-9]' > $temp_out
else
ls -r $base_catalog*$snap_num*0000*particles >> $temp_part
ls -r $base_catalog*$snap_num*0000*halos >> $temp_halo
ls -r $base_catalog*$snap_num*0000*particles | grep -o '_[0-9][0-9][0-9]' >> $temp_out
fi

done

fi

#more $temp_part | tail -n 10

n_files=`wc -l < $temp_part`
echo $n_files

echo 'Starting job on date: ' 
date

# Parallel executable
echo '*******************' 
echo $mpi_exec -n $N_proc $executable $N_read $n_files $output_file $temp_dir $temp_out $temp_part $temp_halo $BOX
echo '*******************' 
$mpi_exec -n $N_proc $executable $N_read $n_files $output_file $temp_dir $temp_out $temp_part $temp_halo $BOX

ls 'output/'$ice_seed'_'$gin_seed'/'*'0000'*'mtree' | wc -l
ls 'output/'$ice_seed'_'$gin_seed'/'*'0000'*'mtree'
