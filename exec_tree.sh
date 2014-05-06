#!/bin/bash

# Merger Tree directories and catalogue path
home_dir='/home/carlesi/'
mtree=$home_dir'/MERGER_TREE/'

# Use serial (1) mpi version (any other number)
use_serial=0

if [ $use_serial -eq 1 ]
then

echo 'Using Serial MTree'

exec_serial=$mtree'/bkp/ahf/bin/MergerTree'
#exec_serial=$mtree'/bkp/ahf/bin/MergerTreeSerial'
#exec_serial=$mtree'/bkp/ahf/bin/MergerTreeSerialNoBoth'
#exec_serial=$mtree'/bkp/MergerTree_v0'

# List two catalogues to test the serial version
#in_1='/home/carlesi/MERGER_TREE/CATALOGUES/SmallSimuMerged/merged_031.AHF_particles'
#in_2='/home/carlesi/MERGER_TREE/CATALOGUES/SmallSimuMerged/merged_030.AHF_particles'
#in_3='/home/carlesi/MERGER_TREE/CATALOGUES/SmallSimuMerged/merged_029.AHF_particles'
in_1='/home/carlesi/MERGER_TREE/CATALOGUES/SussingCatalogs/62.5_dm_061.z0.000.AHF_particles'
in_2='/home/carlesi/MERGER_TREE/CATALOGUES/SussingCatalogs/62.5_dm_060.z0.020.AHF_particles'
in_3='/home/carlesi/MERGER_TREE/CATALOGUES/SussingCatalogs/62.5_dm_059.z0.041.AHF_particles'
out=$mtree'output/serial_sussing_'
#out=$mtree'output/serial_test_'

# Execute the serial version of MergerTree
echo $in_1 
echo $in_2 
echo $in_3 
echo $out

$exec_serial

else 

# N_proc sets the number of mpi tasks
N_proc=1

# Sim Type (0 = test multi snap, 1 = CLUES, 2 = SussingMT 2013)
sim_type=2

# AHF settings
ahf=$home_dir'/MERGER_TREE/latest_ahf/'
executable=$ahf'/bin/MergerTree'

# Store temporary files
temp_dir=$mtree'/temp/'

# Where to dump the lists of the _halo and _particles files
temp_part=$mtree'/temp/list_part_files.ahf'
temp_halo=$mtree'/temp/list_halo_files.ahf'
temp_out=$mtree'/temp/list_files_numbers.ahf'

# If using multiple catalogues (local settings)
if [ $sim_type -eq 0 ]
then

# Box size in mpc (needed for the haloes at the box boundaries, for periodic boxes)
BOX=50000

# Chunks of file
N_read=4

# Output files directory + prefix
output=$mtree'/output/mtree_mpi'
base_catalog=$mtree'/CATALOGUES/SmallSimu/snapshot_'
ls -r $base_catalog*0000*particles > $temp_part
ls -r $base_catalog*0000*halos > $temp_halo
ls -r $base_catalog*0000*particles | grep -o '_0[0-9][0-9]' > $temp_out
fi

if [ $sim_type -eq 1 ]
then

# Box size in kpc (needed for the haloes at the box boundaries, for periodic boxes)
BOX=64000

# Chunks of file
N_read=1

# Output files directory + prefix
output=$mtree'/output/mtree_clues'
base_catalog=$mtree'/clues/4096'
ls -r $base_catalog*particles > $temp_part
ls -r $base_catalog*halos > $temp_halo
ls -r $base_catalog*particles | grep -o '_[0-9][0-9][0-9]' > $temp_out
fi

if [ $sim_type -eq 2 ]
then

# Box size in kpc (needed for the haloes at the box boundaries, for periodic boxes)
BOX=62500

# Chunks of file
N_read=1

# Output files directory + prefix
output=$mtree'/output/mtree_s2013_onetask_dist3mpc'
base_catalog=$mtree'/CATALOGUES/SussingCatalogs/62.5_dm'
ls -r $base_catalog*particles > $temp_part
ls -r $base_catalog*halos > $temp_halo
ls -r $base_catalog*particles | grep -o '_[0-9][0-9][0-9]' > $temp_out
fi


# Check the temporary files, n_files sets the number of files per catalogue to be read in by a single task
#echo $temp_halo
#more $temp_halo
#echo $temp_part
#more $temp_part

n_files=`wc -l < $temp_part`
echo $n_files

cd $ahf
make MergerTree;

# Parallel executable
echo '*******************' 
echo mpiexec -n $N_proc $executable $N_read $n_files $output $temp_dir $temp_out $temp_part $temp_halo $BOX
echo '*******************' 
mpiexec.openmpi -n $N_proc $executable $N_read $n_files $output $temp_dir $temp_out $temp_part $temp_halo $BOX
fi
