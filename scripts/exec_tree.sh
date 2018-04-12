#!/bin/bash

# Merger Tree directories and catalogue path
home_dir='/home/uam37/uam37526/scratch/CLUES/'
mtree=$home_dir'/merger_tree/'
mpi_exec='mpiexec'

# N_proc sets the number of mpi tasks
N_proc=$1

# Chunks of file
N_read=$2

echo $N_proc $N_read

# AHF settings
executable=$mtree'/MergerTree'

# Store temporary files
temp_dir=$mtree'/temp/'

# Where to dump the lists of the _halo and _particles files
temp_part=$mtree'/temp/list_part_files.ahf'
temp_halo=$mtree'/temp/list_halo_files.ahf'
temp_out=$mtree'/temp/list_files_numbers.ahf'

# Box size in kpc (needed for the haloes at the box boundaries, for periodic boxes)
BOX=100000
mkdir $mtree'/output/512_'$gin_seed'/'

run_start=55
run_end=58

for (( irun=$run_start; irun<$run_end; irun++ ))
do

run=`printf %02d ${irun}`

# Output files directory + prefix
mkdir $mtree'/output/512_'$gin_seed'/'$run'/'

output_dir=$mtree'/output/512_'$gin_seed'/'$run'/'
#base_catalog=$mtree'/files/512_'$gin_seed'/'$run'/'
base_catalog=$mtree'/../ahf_catalogs/'$gin_seed'/'$run'/'

output_file=$mtree'/output/512_'$gin_seed'/'$run'/bothways'

if [ "$2" -eq "1" ] ; then 

ls -r $base_catalog*particles > $temp_part
ls -r $base_catalog*halos > $temp_halo
ls -r $base_catalog*particles | grep -o '_[0-9][0-9][0-9]' > $temp_out

else 

ls -r $base_catalog*0000*particles > $temp_part
ls -r $base_catalog*0000*halos > $temp_halo
ls -r $base_catalog*0000*particles | grep -o '_[0-9][0-9][0-9]' > $temp_out

fi

more $temp_part 
# | tail -n 10

n_files=`wc -l < $temp_part`
echo $n_files

echo 'Starting job on date: ' 
date

# Parallel executable
echo '*******************' 
echo $mpi_exec -n $N_proc $executable $N_read $n_files $output_file $temp_dir $temp_out $temp_part $temp_halo $BOX
echo '*******************' 
$mpi_exec -n $N_proc $executable $N_read $n_files $output_file $temp_dir $temp_out $temp_part $temp_halo $BOX
done
