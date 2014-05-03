#!/bin/bash

# N_proc sets the number of mpi tasks
N_proc=4
N_read=1

# AHF settings
home_dir='/home/carlesi/'
ahf=$home_dir'/MERGER_TREE/latest_ahf/'
executable=$ahf'/bin/MergerTree'
exec_serial=$home_dir'/AHF/ahf-v1.0-071/bin/MergerTree'

# Merger Tree directories and catalogue path
mtree=$home_dir'/MERGER_TREE/'
output=$mtree'/output/mtree'
temp_dir=$mtree'/temp/'
temp_part=$mtree'/temp/list_part_files.ahf'
temp_halo=$mtree'/temp/list_halo_files.ahf'

temp_out=$mtree'/temp/list_files_numbers.ahf'
base_out=$mtree'/clues/4096'

# List two catalogues to test the serial version
in_1='/home/carlesi/MERGER_TREE/CATALOGUES/merged_030.AHF_particles'
in_2='/home/carlesi/MERGER_TREE/CATALOGUES/merged_031.AHF_particles'
base_in='/home/carlesi/MERGER_TREE/CATALOGUES/64/merged_'
out=$output'/test_full_30-31.mtr'

#ls -r $base_out*0000*particles > $temp_part
#ls -r $base_out*0000*halos > $temp_halo
#ls -r $base_out*0000*particles | grep -o '_0[0-9][0-9]' > $temp_out

ls -r $base_out*particles > $temp_part
ls -r $base_out*halos > $temp_halo
ls -r $base_out*particles | grep -o '_0[0-9][0-9]' > $temp_out

#echo cd $ahf
cd $ahf
#make clean
make MergerTree;

# Check the temporary files, n_files sets the number of files per catalogue to be read in by a single task
echo $temp_halo
more $temp_halo
echo $temp_part
more $temp_part
n_files=`wc -l < $temp_part`
echo $n_files

# Execute the serial version of MergerTree
#echo $in_1 
#echo $in_2 
#echo $out
#ls $base_in*
#$exec_serial

# Parallel executable
echo '*******************' 
echo mpiexec -n $N_proc $executable $N_read $n_files $output $temp_dir $temp_out $temp_part $temp_halo
echo '*******************' 
mpiexec.openmpi -n $N_proc $executable $N_read $n_files $output $temp_dir $temp_out $temp_part $temp_halo
