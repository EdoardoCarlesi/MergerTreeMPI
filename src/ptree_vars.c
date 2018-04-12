#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "param.h"
#include "tdef.h"

#include "ptree_io.h"
#include "ptree_vars.h"
#include "ptree_funcs.h"
#include "param.h"


#ifdef MTREE_BOTH_WAYS
MTREEptr  *mtree_tmp;	
extern uint64_t  *ncroco_simu;	
#endif

double		Mparticle;	
int		totComp=0;

uint64_t    nHalos[3], nHalosTmp[2];

size_t      totHaloSize, totHaloSizeTmp;	

double BoxSize; 	

int NReadTask, TotTask, LocTask; 
int SendTask, RecvTask, DeltaTask;
int filesPerTask, extraFilesPerTask;
MPI_Status status;

HALOptr halos[3];
HALOptr halos_tmp[3]; 		
HALO_MPIptr halos_mpi;

/* Buffers to recieve and swap the new halos across the tasks */
uint64_t nHalosBuffer[3];	
uint64_t nHalosSendBuffer[3];	
uint64_t nHalosRecvBuffer[3];	
MTREEptr  *mtree_tmp;	
uint64_t  *ncroco_simu;	




