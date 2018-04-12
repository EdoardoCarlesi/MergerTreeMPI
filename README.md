MergerTreeMPI
=============
Edoardo Carlesi - 2013
edoardo.carlesi@uam.es


This is a massively MPI parallel version of the original MergerTree.c code, part of the AHF analysis package.
The core algorithm of the old serial code has been substantially changed to be able to run on shared memory machines.
In fact,the identification of the particles' position within an array with their ID is no longer possible, as each
task now stores a different array where the particles ID do not necessarily match their positions within it.
This reduces substantially the memory requirements while on the other hand slows down the process of comparison, which
however scales as Ntask^2.
