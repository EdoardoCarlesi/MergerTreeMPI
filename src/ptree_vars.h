#include <mpi.h>
#include <stdlib.h>

#include "param.h"
#include "tdef.h"
//#include "libutility/utility.h"

#ifndef PTREE_VARS_H
#define PTREE_VARS_H

typedef struct MTREE 
{
  uint64_t haloid[2];
  uint64_t id[2];
  uint64_t npart[2];
  uint64_t common;
} MTREE;

typedef struct MTREE *MTREEptr;

typedef struct HALOS
{
  uint64_t  haloid;
  uint64_t  npart;
  float    Xc[3];		/* C.o.m. position of the halo */
  float    Vc[3];		/* C.o.ma. position of the halo */
  uint64_t *Pid;		/* Pid holds the global particle ID */
  uint64_t  ncroco;		/* Local number of cross correlation */
  uint64_t  global_ncroco;	/* Global number - sums all the haloes swapped so far */
  MTREEptr  mtree;
} HALOS;

typedef struct HALOS *HALOptr;

/* Structure used for MPI communication */
typedef struct HALO_MPI *HALO_MPIptr;

typedef struct HALO_MPI
{
  uint64_t  haloid;
  uint64_t  npart;
  uint64_t  ncroco;			/* Global number - sums all the haloes swapped so far */
  uint64_t  haloid_max_merit;		/* ID of the descendant halo with max_merit */
  uint64_t  id_max_merit;		/* Local HALO position of the descendant halo with max_merit */
} HALO_MPI;

extern HALOptr halos[3];
extern HALOptr halos_tmp[3]; 		
extern HALO_MPIptr halos_mpi;

/* Buffers to recieve and swap the new halos across the tasks */
extern uint64_t nHalosBuffer[3];	
extern uint64_t nHalosSendBuffer[3];	
extern uint64_t nHalosRecvBuffer[3];	

//#ifdef MTREE_BOTH_WAYS
extern MTREEptr  *mtree_tmp;	/* When doing the mtree both ways, we need to store the mtree for the 1->0, since
				   simu 1 is the one which is being swapped through the tasks */
extern uint64_t  *ncroco_simu;	/* Keep track of the total number of croco in simu1 on the local task */
//#endif

extern double		Mparticle;	
extern int		totComp;

extern uint64_t    nHalos[3];
extern uint64_t    nHalosTmp[2];

/* When reading-in the halo keep track of its size including dynamically allocated memory */
extern size_t      totHaloSize;	
extern size_t      totHaloSizeTmp;	

/* Needed for periodic boundary conditions */
extern double BoxSize; 	

/*	MPI Communication variables	*/
extern int NReadTask;		
extern int TotTask; 
extern int LocTask; 
extern int SendTask; 
extern int RecvTask;
extern int DeltaTask;
extern int filesPerTask; 
extern int extraFilesPerTask;
extern MPI_Status status;
#endif
