#include <mpi.h>
#include "ptree_vars.h"

#ifndef PTREE_FUNCS_H
#define PTREE_FUNCS_H

int	cross_correlation      (int  iloop, int simu0, int simu1); 
int	create_mtree           (uint64_t ihalo, int isimu0, int isimu1, int iloop);
int	clean_connection       (uint64_t ihalo, int isimu0, int isimu1, int iloop);
int	order_by_merit		(int isimu0);
uint64_t max_merit              (uint64_t ihalo);

int	bubblesort			(uint64_t *array, uint64_t size);
int	selectionsort			(uint64_t *array, uint64_t size);
uint64_t	cmpfunc			(uint64_t *a, uint64_t *b); 
uint64_t	cmpdbl				(double *a, double *b); 

int	add_halos			(int ifile, int isimu);
int	copy_halos			(int isimu0, int isimu1);

int	copy_halos_to_halos_mpi		(int isimu0);
int	MPI_Swaphalos			(int isimu);
int	MPI_Swaphalos_mpi		(int isimu);

int	alloc_halos			(int isimu);
int	free_halos			(int isimu);
int	order_halo_ids			(int isimu);

int	load_balance			(int isimu); 
void	intersection 			(int isimu0, int isimu1, uint64_t ihalo, uint64_t khalo, uint64_t *common);

int	compute_com_distance		(uint64_t ihalo0, uint64_t ihalo1, int isimu0, int isimu1);

uint64_t	constructHaloId		(uint64_t ihalo);

void 	check_halos			(int isimu);
#endif
