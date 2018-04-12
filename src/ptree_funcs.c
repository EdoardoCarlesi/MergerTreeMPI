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

#include "libutility/utility.h"


/*==================================================================================================
 * load_balance:
 * 
 * when using more tasks than input files (parts of a given catalogue) we need to redistribute
 * halos among them. 
 * At this point only the halos structures have been read-in; the inverse
 * mapping to the particles will be done after the halos have been scattered.
 * Since tasks that read in the halos use dynamically allocated memory (for Pid, that is
 * relative to the particles belonging to the halo) we use the MPI_Pack/MPI_Unpack functions.
 *
 *==================================================================================================*/
int load_balance(int isimu)
{
  int irecv=0, itask=0, nTaskPerSend=0, extraTask=0, buffer_pos_recv=0;
  uint64_t ihalo=0, jhalo=0;
  int *nRecvTasks=NULL, *sendTasks=NULL, **recvTasks=NULL, *nHalosBuffer=NULL;
  uint64_t locNPart=0; 
  void **haloSendBuffer=NULL, *haloRecvBuffer=NULL;

  size_t sizeUInt64=0, sizeRecvBuffer=0, sizeTempBuffer=0, sizeLocHalo=0; 
  size_t *sizeSendBuffer=NULL; 

  HALOptr temp_halo[2];
  int *buffer_position=NULL;

#ifdef IMPROVE_LB
  uint64_t *partPerTask=NULL, minNPartTmp=0, minNPart=0;
  double nHaloFac=0, nPartFac=0, nPartFacTmp=0;
#endif

  sizeUInt64 = sizeof(uint64_t);

  /* recv tasks keeps track of all the tasks to which a single (reading) tasks 
   * has to broadcast data, send task tells the recieving task who is sending.
   * Each reading tasks can send to many tasks, but each recieving task
   * only gets data from a single one. */
  if(LocTask < NReadTask)
  { 
    recvTasks = (int **) calloc(NReadTask, sizeof(int *));
    nRecvTasks = (int *) calloc(NReadTask, sizeof(int));
  }
  else
  {
    sendTasks = (int *) calloc(DeltaTask, sizeof(int));
  }

  nTaskPerSend = (int) DeltaTask / NReadTask;
  extraTask = DeltaTask % NReadTask;

  if(LocTask == 0)
    fprintf(stderr, "Load balancing from %d ReadTasks on %d TotalTasks, each SendTask has %d RecvTasks\n.",
             NReadTask, TotTask, nTaskPerSend);

  /* we figure out who is sending to how many tasks */ 
  if(LocTask < NReadTask)
  {
     if(LocTask < extraTask)
     {
	nRecvTasks[LocTask] = nTaskPerSend + 1;
        recvTasks[LocTask] = (int *) calloc(nTaskPerSend+1, sizeof(int));
     }
     else if (LocTask >= extraTask && nTaskPerSend > 0)
     {
	nRecvTasks[LocTask] = nTaskPerSend;
        recvTasks[LocTask] = (int *) calloc(nTaskPerSend, sizeof(int));
     }
     else 
     {
	nRecvTasks[LocTask] = 0;
        recvTasks[LocTask] = NULL;
     }

     for(itask=0; itask<nRecvTasks[LocTask]; itask++)
     {
       recvTasks[LocTask][itask] = LocTask + (itask + 1) * NReadTask;
     }
   } 
   else /* the recieving tasks have to figure out who is sending */ 
   {
      sendTasks[LocTask-NReadTask] = LocTask % NReadTask;
   }

  /* Figure out how much every task has to send */

  if(LocTask < NReadTask && nRecvTasks[LocTask]>0)
  {
    haloSendBuffer = (void **) calloc(nRecvTasks[LocTask], sizeof(void *)); 
    sizeSendBuffer = (size_t *) calloc(nRecvTasks[LocTask], sizeof(size_t));
    buffer_position = (int *) calloc(nRecvTasks[LocTask], sizeof(int));
    nHalosBuffer = (int *) calloc(nRecvTasks[LocTask]+1, sizeof(int));

    haloSendBuffer[irecv] = (void *) malloc(1);
 
    temp_halo[isimu] = (HALOptr) calloc(1, sizeof(HALOS)); 

#ifdef IMPROVE_LB
    partPerTask = (uint64_t *) calloc(nRecvTasks[LocTask]+1, sizeUInt64); 
#endif

    /* Loop over haloes and assign them to the different tasks, packing them into separate buffers */
    for(ihalo=1; ihalo<=nHalos[isimu]; ihalo++)
    {
#ifdef IMPROVE_LB
      jhalo = ihalo-1;
#else
      jhalo = nHalos[isimu] - ihalo;
#endif

      irecv = (jhalo) % (nRecvTasks[LocTask]+1);
      locNPart = halos[isimu][jhalo].npart;
      sizeLocHalo = (size_t) locNPart * sizeUInt64;
     
#ifdef IMPROVE_LB
      minNPart = partPerTask[0];
      nHaloFac = sqrt(nHalosBuffer[0]);
      nPartFac = nHaloFac * (double) minNPart;
      
      for(itask=1; itask<nRecvTasks[LocTask]+1; itask++)
      {
         minNPartTmp = partPerTask[itask];
 	 nPartFacTmp = (double) minNPartTmp * LB_PART_FAC * sqrt(nHalosBuffer[itask] * LB_HALO_FAC);

         if(nPartFacTmp < nPartFac) 
         {
	    irecv = itask;
	    nPartFac = nPartFacTmp;
         }
      }

         partPerTask[irecv] += locNPart;

#endif

      /* Each halo holds two uint64_t, three Xc floats and two arrays Pid */
      sizeTempBuffer = sizeUInt64 * (halos[isimu][jhalo].npart + 2) + 3 * sizeof(float);

      if(irecv == nRecvTasks[LocTask]) /* keep the halo locally */
      {
         temp_halo[isimu] = (HALOptr) realloc(temp_halo[isimu], (nHalosBuffer[irecv] + 1) * sizeof(HALOS));
         temp_halo[isimu][nHalosBuffer[irecv]].npart = locNPart;
         temp_halo[isimu][nHalosBuffer[irecv]].haloid = halos[isimu][jhalo].haloid;

         temp_halo[isimu][nHalosBuffer[irecv]].Xc[0] = halos[isimu][jhalo].Xc[0];
         temp_halo[isimu][nHalosBuffer[irecv]].Xc[1] = halos[isimu][jhalo].Xc[1];
         temp_halo[isimu][nHalosBuffer[irecv]].Xc[2] = halos[isimu][jhalo].Xc[2];

         temp_halo[isimu][nHalosBuffer[irecv]].Pid = (uint64_t *) malloc(sizeLocHalo);

         memcpy(temp_halo[isimu][nHalosBuffer[irecv]].Pid, halos[isimu][jhalo].Pid, sizeLocHalo);   
      }
      else
      {
         /* now allocate and initialize the buffers to which will be sent to the next task */
         sizeSendBuffer[irecv] += sizeTempBuffer;
         haloSendBuffer[irecv]  = (uint64_t *) realloc(haloSendBuffer[irecv], sizeSendBuffer[irecv]);

         /* Now pack the halo and particles information into a single buffer */
         MPI_Pack(&halos[isimu][jhalo].haloid, sizeUInt64, MPI_BYTE, haloSendBuffer[irecv],
          sizeSendBuffer[irecv], &buffer_position[irecv], MPI_COMM_WORLD);

         MPI_Pack(&halos[isimu][jhalo].npart, sizeUInt64, MPI_BYTE, haloSendBuffer[irecv],
          sizeSendBuffer[irecv], &buffer_position[irecv], MPI_COMM_WORLD);

         MPI_Pack(&halos[isimu][jhalo].Xc[0], 3 * sizeof(float), MPI_BYTE, haloSendBuffer[irecv],
          sizeSendBuffer[irecv], &buffer_position[irecv], MPI_COMM_WORLD);

         MPI_Pack(&halos[isimu][jhalo].Pid[0], sizeLocHalo, MPI_BYTE, haloSendBuffer[irecv], 
          sizeSendBuffer[irecv], &buffer_position[irecv], MPI_COMM_WORLD);
      }

      nHalosBuffer[irecv]++;

      /* reallocate the local halos struct at each step to free up some memory */
#ifndef IMPROVE_LB
      free(halos[isimu][jhalo].Pid);
      halos[isimu] = (HALOptr) realloc(halos[isimu], (jhalo+1) * sizeof(HALOS));
#endif
   } /* for ihalo */

   nHalos[isimu] = nHalosBuffer[nRecvTasks[LocTask]];

   free(halos[isimu]);
   halos[isimu] = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));

   for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
   { 
      halos[isimu][ihalo].npart = temp_halo[isimu][ihalo].npart;
      halos[isimu][ihalo].haloid = temp_halo[isimu][ihalo].haloid;

      halos[isimu][ihalo].Xc[0] = temp_halo[isimu][ihalo].Xc[0];
      halos[isimu][ihalo].Xc[1] = temp_halo[isimu][ihalo].Xc[1];
      halos[isimu][ihalo].Xc[2] = temp_halo[isimu][ihalo].Xc[2];

      locNPart = halos[isimu][ihalo].npart;
      sizeLocHalo = (size_t) locNPart * sizeUInt64;

      halos[isimu][ihalo].Pid = (uint64_t *) malloc (sizeLocHalo);

      memcpy(halos[isimu][ihalo].Pid, temp_halo[isimu][ihalo].Pid, sizeLocHalo);   
      
      free(temp_halo[isimu][ihalo].Pid);
   }

#ifdef DEBUG_MPI
   fprintf(stderr, "Task=%d is left with %"PRIu64" haloes.\n", LocTask, nHalos[isimu]);
#endif

   free(temp_halo[isimu]);

#ifdef DEBUG_LOG 
   if(LocTask < NReadTask)
   {
      for(khalo=0; khalo<nHalos[isimu]; khalo++)
	dump_log_halo(isimu, khalo);
   }
#endif

     if(nHalos[isimu] == 0)
     {
       fprintf(stderr, "WARNING! No haloes left on Task=%d. Are you using too many MPI tasks?\n", LocTask);
     }
     else
     {
#ifdef DEBUG_MPI
	for(itask=0; itask<nRecvTasks[LocTask]; itask++)
 	{
	    fprintf(stderr, "Task=%d is sending data to task %d\n",
	            LocTask, recvTasks[LocTask][itask]);
	}

	    fprintf(stderr, "Task=%d has done packing data, local data has %"PRIu64" halos.\n", 
		    LocTask, nHalos[isimu]);
#endif
     }
  } /* if(LocTask < NReadTask)*/

  /* Send all the buffers */
  if(LocTask < NReadTask && nRecvTasks[LocTask]>0)
  {
    for(irecv=0; irecv<nRecvTasks[LocTask]; irecv++)
    {
      MPI_Send(&sizeSendBuffer[irecv], sizeof(size_t), MPI_BYTE, recvTasks[LocTask][irecv], 0, MPI_COMM_WORLD);
      MPI_Send(&nHalosBuffer[irecv], 1, MPI_INT, recvTasks[LocTask][irecv], 0, MPI_COMM_WORLD);
      MPI_Send(haloSendBuffer[irecv], sizeSendBuffer[irecv], MPI_BYTE, recvTasks[LocTask][irecv], 0, MPI_COMM_WORLD);

#ifdef DEBUG_MPI
    fprintf(stderr, "Sending messages from task=%d to task=%d.\n", LocTask, recvTasks[LocTask][irecv]);
#endif

    }
  }
  else if(LocTask >= NReadTask) /* the task is a recieving one */
  {
      MPI_Recv(&sizeRecvBuffer, sizeof(size_t), MPI_BYTE, sendTasks[LocTask-NReadTask], 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&nHalos[isimu], 1, MPI_INT, sendTasks[LocTask-NReadTask], 0, MPI_COMM_WORLD, &status);
   
      /* Use this buffer to store the incoming packed halos */
      haloRecvBuffer = (void *) malloc(sizeRecvBuffer);

#ifdef DEBUG_MPI
    fprintf(stderr, "\nRecieved %"PRIu64" halos on task=%d from task=%d.\n", nHalos[isimu], LocTask, sendTasks[LocTask-NReadTask]);
#endif

      MPI_Recv(haloRecvBuffer, sizeRecvBuffer, MPI_BYTE, sendTasks[LocTask-NReadTask], 0, MPI_COMM_WORLD, &status);
  }

  /* Unpack the recieved buffers on each task */
  if(LocTask >= NReadTask)
  {
    halos[isimu]   = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));
    buffer_pos_recv = 0;

    for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
    {
       /* Now unpack halos particles information from the buffers into the halos structures */
       MPI_Unpack(haloRecvBuffer, sizeRecvBuffer, &buffer_pos_recv, &halos[isimu][ihalo].haloid, 
	sizeUInt64, MPI_BYTE, MPI_COMM_WORLD);
	
       MPI_Unpack(haloRecvBuffer, sizeRecvBuffer, &buffer_pos_recv, &halos[isimu][ihalo].npart, 
 	sizeUInt64, MPI_BYTE, MPI_COMM_WORLD);

       MPI_Unpack(haloRecvBuffer, sizeRecvBuffer, &buffer_pos_recv, &halos[isimu][ihalo].Xc[0], 
 	3 * sizeof(float), MPI_BYTE, MPI_COMM_WORLD);

       locNPart = halos[isimu][ihalo].npart;
       sizeLocHalo = (size_t) locNPart * sizeUInt64;

       halos[isimu][ihalo].Pid = (uint64_t *) malloc(sizeLocHalo);

       MPI_Unpack(haloRecvBuffer, sizeRecvBuffer, &buffer_pos_recv, &halos[isimu][ihalo].Pid[0], 
        sizeLocHalo, MPI_BYTE, MPI_COMM_WORLD);
    }

   } /* if(LocTask >= NReadTask)*/

#ifdef DEBUG_MPI
	size_t totHaloSizeTest=0;

	for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
	  totHaloSizeTest += 2 * (1 + halos[isimu][ihalo].npart) * sizeUInt64;

	fprintf(stderr, "\nDone load balancing. On Task=%d there are %"PRIu64" halos, total size=%zd MB.\n", 
		LocTask, nHalos[isimu], totHaloSizeTest/1024/1024);
#endif

  /* Now free local buffers */
  if(LocTask < NReadTask && nRecvTasks[LocTask]>0)
  {
     for(ihalo=0; ihalo<nRecvTasks[LocTask]; ihalo++)
       free(haloSendBuffer[ihalo]);
     
     free(haloSendBuffer);
     free(sizeSendBuffer);
     free(nHalosBuffer);
  } 
  else
  {
     free(haloRecvBuffer);
  }
	
  return(1);
}


/* This function sorts the Pids in ascending order, to enable a faster comparison of the halo contents */
int order_halo_ids(int isimu)
{
   uint64_t ihalo=0, locNPart=0, *Pidord=NULL, ipart=0, iord=0;
   double *Pidnew=NULL;
 
   /* Store the ordered IDs for successive comparison */
   for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
   {

      locNPart = halos[isimu][ihalo].npart;
      Pidord = calloc(locNPart, sizeof(uint64_t));
      memcpy(Pidord, halos[isimu][ihalo].Pid, locNPart * sizeof(uint64_t));   

      qsort(Pidord, locNPart, sizeof(uint64_t), &cmpfunc);
	
	iord = 0;
	for (ipart=1; ipart<locNPart; ipart++)
	{
		if (Pidord[0] > Pidord[ipart])
			iord++;
	}
		
		// This check is necessary as qsort is not stable when using very long uint64 ids
		if (iord > 0 && ihalo < 10)
			fprintf(stderr, "Error in qsort for halo %"PRIu64" npart: %"PRIu64" in simu: %"PRIu64" ord:%"PRIu64", resorting IDs.\n", 
				ihalo, locNPart, isimu, iord);

		if (iord > 0)
		{	
			// The ids are converted to double, this allows a more precise comparison.
      			Pidnew = calloc(locNPart, sizeof(double));
			
			for (ipart=0; ipart<locNPart; ipart++)
			{
				Pidnew[ipart] = (double) Pidord[ipart];
				Pidnew[ipart] /= 1000.0;
			}

      			qsort(Pidnew, locNPart, sizeof(double), &cmpdbl);

			for (ipart=0; ipart<locNPart; ipart++)
			{
				Pidnew[ipart] *= 1000.0;
				Pidord[ipart] = (uint64_t) Pidnew[ipart];
			}

			iord=0;
			for (ipart=0; ipart<locNPart; ipart++)
			{
				if (Pidord[0] > Pidord[ipart])
				iord++;
			}


			// SANITY CHECK 
			if (iord > 0 && ihalo < 10)
			{
				fprintf(stderr, "Again something went wrong for halo %"PRIu64" npart: %"PRIu64" in simu: %"PRIu64" ord:%"PRIu64".\n", 
				ihalo, locNPart, isimu, iord);
					
				for (ipart=0; ipart<10; ipart++)
					fprintf(stderr, "ord ID[%"PRIu64"] = %"PRIu64", %lf \n", ipart, Pidord[ipart], Pidnew[ipart]);

			}

       			if (Pidnew != NULL)
       			{
        			 free(Pidnew);
        			 Pidnew=NULL;
       			}
		}

      halos[isimu][ihalo].Pid = Pidord;
    }
 
  return 1;
}


/*
 * This function swaps the halo structures across the tasks using the MPI_Pack / MPI_Unpack functions 
 * to be able to send dynamically allocated data.
 */
int MPI_Swaphalos(int isimu)
{
  int buffer_pos_recv=0, buffer_position=0;
  uint64_t ihalo=0; 
  size_t sizeUInt64=0, sizeSendBuffer=0, sizeRecvBuffer=0, sizeLocHalo=0; 
  void *haloSendBuffer=NULL, *haloRecvBuffer=NULL;
  uint64_t locNPart=0; 

  sizeUInt64 = sizeof(uint64_t);

#ifdef DEBUG_MPI
    fprintf(stderr, "MPI_Swapping halo[%d] to task=%d from task=%d.\n", isimu, LocTask, SendTask);
#endif

    /* Init buffers */
//    haloSendBuffer = (void *) malloc(1); 
//    haloRecvBuffer = (void *) malloc(1); 

    /* Loop over haloes and assign them to the different tasks, packing them into separate buffers */
    for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
    {
      locNPart = halos[isimu][ihalo].npart;
      sizeLocHalo = locNPart * sizeUInt64; 

      /* Each halo holds two uint64_t, and array Pid and the coordinates */
      sizeSendBuffer += sizeUInt64 * (locNPart + 2) + 3 * sizeof(float);

	if(ihalo == 0)  haloSendBuffer = (void *) malloc(sizeSendBuffer);
      	else 		haloSendBuffer = (void *) realloc(haloSendBuffer, sizeSendBuffer);

      /* Pack the halo and particles information into a single buffer */
      MPI_Pack(&halos[isimu][ihalo].haloid, sizeUInt64, MPI_BYTE, haloSendBuffer,
        sizeSendBuffer, &buffer_position, MPI_COMM_WORLD);

      MPI_Pack(&halos[isimu][ihalo].npart, sizeUInt64, MPI_BYTE, haloSendBuffer,
        sizeSendBuffer, &buffer_position, MPI_COMM_WORLD);

      MPI_Pack(&halos[isimu][ihalo].Xc[0], 3 * sizeof(float), MPI_BYTE, haloSendBuffer,
        sizeSendBuffer, &buffer_position, MPI_COMM_WORLD);

      MPI_Pack(&halos[isimu][ihalo].Pid[0], sizeLocHalo, MPI_BYTE, haloSendBuffer, 
        sizeSendBuffer, &buffer_position, MPI_COMM_WORLD);
	
      /* Free Pid space once it's packed */
      free(halos[isimu][ihalo].Pid);
	/* At some later stages the mtree might be allocated already */
      if(halos[isimu][ihalo].mtree != NULL)
	      free(halos[isimu][ihalo].mtree);
    } /* for ihalo */

   free(halos[isimu]);

#ifdef DEBUG_MPI
   fprintf(stderr, "nhalos to be sent %"PRIu64" loc_send_size = %zd bytes, loc=%d, recv=%d send=%d\n", 
		nHalos[isimu], sizeSendBuffer, LocTask, RecvTask, SendTask);
#endif

   MPI_Sendrecv(&nHalos[isimu], sizeof(uint64_t), MPI_BYTE, RecvTask, 0,
                &nHalosBuffer[isimu],  sizeof(uint64_t), MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

   MPI_Sendrecv(&sizeSendBuffer, sizeof(size_t), MPI_BYTE, RecvTask, 0,
                &sizeRecvBuffer, sizeof(size_t), MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

#ifdef DEBUG_MPI
   fprintf(stderr, "nhalos to recieve %"PRIu64" loc_recv_size = %zd bytes, loc=%d, recv=%d send=%d\n", 
		nHalosBuffer[isimu], sizeRecvBuffer, LocTask, RecvTask, SendTask);
#endif

   haloRecvBuffer = (void *) malloc(sizeRecvBuffer);
   nHalos[isimu] = nHalosBuffer[isimu];

   MPI_Sendrecv(haloSendBuffer, sizeSendBuffer, MPI_BYTE, RecvTask, 0,
                haloRecvBuffer, sizeRecvBuffer, MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

    free(haloSendBuffer);

#ifdef DEBUG_MPI
    fprintf(stderr, "Task=%d has done packing data, local data size is now=%zd, with %"PRIu64" halos.\n", 
	LocTask, sizeSendBuffer, nHalos[isimu]);

    fprintf(stderr, "Sending %zd bytes of messages on task=%d from task=%d.\n", 
	sizeRecvBuffer, LocTask, SendTask);
#endif

    /* Now unpack halos particles information from the buffers into the halos structures */
    halos[isimu] = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));
    buffer_pos_recv = 0;

    for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
    {
       MPI_Unpack(haloRecvBuffer, sizeRecvBuffer, &buffer_pos_recv, &halos[isimu][ihalo].haloid, 
	sizeUInt64, MPI_BYTE, MPI_COMM_WORLD);
	
       MPI_Unpack(haloRecvBuffer, sizeRecvBuffer, &buffer_pos_recv, &halos[isimu][ihalo].npart, 
 	sizeUInt64, MPI_BYTE, MPI_COMM_WORLD);

       MPI_Unpack(haloRecvBuffer, sizeRecvBuffer, &buffer_pos_recv, &halos[isimu][ihalo].Xc[0], 
 	3 * sizeof(float), MPI_BYTE, MPI_COMM_WORLD);

       locNPart = halos[isimu][ihalo].npart;
       sizeLocHalo = (size_t) locNPart * sizeUInt64;

       halos[isimu][ihalo].Pid = (uint64_t *) malloc(sizeLocHalo);
       halos[isimu][ihalo].mtree = NULL; //(MTREEptr) malloc(sizeof(MTREE));

       MPI_Unpack(haloRecvBuffer, sizeRecvBuffer, &buffer_pos_recv, &halos[isimu][ihalo].Pid[0], 
        sizeLocHalo, MPI_BYTE, MPI_COMM_WORLD);
   } /* if(LocTask >= NReadTask)*/

#ifdef DEBUG_MPI
    fprintf(stderr, "Task=%d has recived data, local data size is now=%zd, with %"PRIu64" halos.\n", 
	LocTask, sizeRecvBuffer, nHalos[isimu]);
#endif

  /* Now free local buffer */
  free(haloRecvBuffer); 

  return(1);
}


/*==================================================================================================
 * cross_correlation:
 *
 *  for each halo at isimu=0 figure out how many particles are in common with khalo at isimu=1
 *
 *==================================================================================================*/
int cross_correlation(int iloop, int isimu0, int isimu1)
{
  uint64_t  ihalo, count=COUNT_VAR;
  clock_t   elapsed = (time_t)0;
  
  /*---------------------------------------------------------
   * backwards correlation
   *---------------------------------------------------------*/
  if(LocTask == 0 && iloop == 0)
    fprintf(stderr,"  o generating cross-correlation %d->%d for %"PRIu64" haloes on %d tasks...",
		isimu0, isimu1, nHalos[isimu0], TotTask);

#ifdef DEBUG_MPI
   fprintf(stderr, "\ncreating mtree for %"PRIu64" haloes on task %d\n", nHalos[0], LocTask);
#endif 

    elapsed -= time(NULL);

  /* cross-correlation simu0->simu1. When creating the m_tree now we need to 
   * take into account that simu1 is split in TotTask files, and we now are
   * creating the m_tree for the iloop-th one
   *  */
  for(ihalo=0; ihalo<nHalos[isimu0]; ihalo++) 
  {
#ifdef USE_HALO_MIN_PART
    if(halos[isimu0][ihalo].npart > HALO_MIN_PART)
    {
#endif
        create_mtree(ihalo, isimu0, isimu1, iloop);

#ifdef USE_HALO_MIN_PART
    }
#endif

        if(ihalo == count) 
        {
#ifdef DEBUG_MPI	
           fprintf(stderr, "on task=%d, jpart=%"PRIu64"/%"PRIu64" done.\n", LocTask, ihalo, nHalos[0]);
#else
           fprintf(stderr, ".");
#endif
           count += 500;
        }
 }

    elapsed += time(NULL); 
  
  if(LocTask == 0 && iloop == 0)
	  fprintf(stderr,"\n done in %4.2f sec.\n", (double)elapsed);

  return(1);
}


/*==================================================================================================
 * create_mtree
 *==================================================================================================*/
int create_mtree(uint64_t ihalo, int isimu0, int isimu1, int iloop)
{
  uint64_t  khalo, ncroco, jcroco;
  uint64_t  *common=NULL;
  	
  /* common[] records how many particles ihalo(isimu0) has in common with khalo(isimu1) */
  common = (uint64_t *) calloc(nHalos[isimu1], sizeof(uint64_t));

  for(khalo=0; khalo<nHalos[isimu1]; khalo++) 
  {
    /* Only take into account pairs of haloes closer than MAX_HALO_DIST*/
#ifdef MAX_HALO_DIST
    if (compute_com_distance(ihalo, khalo, isimu0, isimu1) == 1)
    {
#else
	if(ihalo==0 && khalo==0 && LocTask ==0) 
	{
	   fprintf(stderr, "\nComputing halo cross correlation with no distance cutoff.\n");
#endif
     	   	intersection(isimu0, isimu1, ihalo, khalo, common);
        }
  }

  /* determine number of credible cross-correlations */
  ncroco = 0;
  for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
    if(common[khalo] > MINCOMMON)
      ncroco++;
  }

  /* Update the number of connections found in total - ncroco only holds the LOCAL connections */
  halos[isimu0][ihalo].global_ncroco += ncroco; 

  /* does not make sense to continue if there are no cross-correlations */
  if(ncroco > 0) 
  {
    /* allocate memory for cross-correlations */
    if (halos[isimu0][ihalo].mtree == NULL) 
	halos[isimu0][ihalo].mtree = (MTREEptr) calloc(halos[isimu0][ihalo].global_ncroco, sizeof(MTREE)); 
    else 
	halos[isimu0][ihalo].mtree = (MTREEptr) realloc(halos[isimu0][ihalo].mtree, (halos[isimu0][ihalo].global_ncroco) * sizeof(MTREE)); 
  
    /* The halo[][].mtree for iloop>0 may already contain some connections to haloes 
     * so we need to initialize the variables */
    jcroco = halos[isimu0][ihalo].global_ncroco - ncroco; 

    /* Store mtree inside halos[][] structure in incorrect order. 
     * The correct order according to merit function will be implemented when all haloes are stored */
    for(khalo=0; khalo<nHalos[isimu1]; khalo++) 
    {
      if(common[khalo] > MINCOMMON)
      {
        halos[isimu0][ihalo].mtree[jcroco].id[0]    = ihalo;
        halos[isimu0][ihalo].mtree[jcroco].haloid[0]= halos[isimu0][ihalo].haloid;
        halos[isimu0][ihalo].mtree[jcroco].npart[0] = halos[isimu0][ihalo].npart;
        halos[isimu0][ihalo].mtree[jcroco].common   = common[khalo];
        halos[isimu0][ihalo].mtree[jcroco].id[1]    = khalo;
        halos[isimu0][ihalo].mtree[jcroco].haloid[1]= halos[isimu1][khalo].haloid;
        halos[isimu0][ihalo].mtree[jcroco].npart[1] = halos[isimu1][khalo].npart;
        jcroco++;
      }
    }  
  } /* if ncroco>0 */
  
    /* free temporary structures */
    if(common) free(common);

  return(1);
}



#ifdef MAX_HALO_DIST
/* Computes the distance of the center of mass of two given halos. If their distance is larger than MAX
 * then the comparison will be skipped. Assumes kpc units.
 */
int compute_com_distance(uint64_t ihalo0, uint64_t ihalo1, int isimu0, int isimu1)
{
  int Dist = 0, icoord = 0;
  float X0[3], X1[3], dX[3], dbox, dx, dSum, dDist;
#ifdef DYN_MAX_DIST
  double dynDist, M1, M2, Mnorm, r1, r2;
#endif

  for(icoord = 0; icoord < 3; icoord++)
  {
    X0[icoord] = halos[isimu0][ihalo0].Xc[icoord];
    X1[icoord] = halos[isimu1][ihalo1].Xc[icoord];
  }
  
  dSum = 0;

  for(icoord = 0; icoord < 3; icoord++)
  {
    dx = X1[icoord] - X0[icoord];
    dx = sqrt(dx * dx);

    if(X0[icoord] > X1[icoord])
      dbox = (BoxSize - X0[icoord]) + X1[icoord];
    else
      dbox = (BoxSize - X1[icoord]) + X0[icoord];

    if(dx < dbox)
      dX[icoord] = dx;
    else
      dX[icoord] = dbox;

    dSum += pow2(dX[icoord]);
  }

  dDist = sqrt(dSum);

#ifdef DYN_MAX_DIST
  M1 = Mparticle * halos[isimu0][ihalo0].npart;
  M2 = Mparticle * halos[isimu1][ihalo1].npart;
  Mnorm = 1.e+14;

  r1 = (M1/Mnorm); //, 1./3.);
  r2 = (M2/Mnorm); //, 1./3.);

  // This is normalized to a radius of MAX_HALO_DIST for a structure of 1e+14 MSun
  dynDist = MAX_HALO_DIST * (r1 + r2) * MAX_HALO_DIST * MAX_HALO_DIST;

	//fprintf(stderr, "M1=%e M2=%e, r1=%f r2=%f, dynDist=%f\n", M1, M2, r1, r2, dynDist);

  if (dynDist > MAX_HALO_DIST) dynDist = MAX_HALO_DIST;		// MAX HALO DIST is the upper limit	

	// Use CUBES to compare instead of sqrt, it takes much less time
  if(dDist * dDist * dDist < dynDist)
    Dist = 1;
#else
  if(dDist < MAX_HALO_DIST)
    Dist = 1;
#endif

  return Dist;
}
#endif

int copy_halos(int isimu0, int isimu1)
{
   uint64_t ihalo, npart;
   size_t sizeHalo;

   for(ihalo = 0; ihalo<nHalos[isimu1]; ihalo++)
   {
      npart = halos[isimu0][ihalo].npart;
      halos[isimu1][ihalo].npart = npart;
      halos[isimu1][ihalo].haloid = halos[isimu0][ihalo].haloid;
      memcpy(halos[isimu1][ihalo].Xc, halos[isimu0][ihalo].Xc, 3 * sizeof(float));
	
      sizeHalo = npart * sizeof(uint64_t);
 
      halos[isimu1][ihalo].Pid = (uint64_t *) malloc(sizeHalo);
      memcpy(halos[isimu1][ihalo].Pid, halos[isimu0][ihalo].Pid, sizeHalo);   
   }

  return 1;
}


/* We copy the "static" part of the HALOS structure into a HALO_MPI type which can be easily handled in MPI communication */
int copy_halos_to_halos_mpi(int isimu1)
{
  uint64_t ihalo;

  for(ihalo = 0; ihalo<nHalos[isimu1]; ihalo++)
  {	
     halos_mpi[ihalo].haloid = halos[isimu1][ihalo].haloid; 
     halos_mpi[ihalo].npart = halos[isimu1][ihalo].npart; 
     halos_mpi[ihalo].ncroco = halos[isimu1][ihalo].global_ncroco; 
	
	if(halos[isimu1][ihalo].mtree == NULL)
	{
		halos_mpi[ihalo].haloid_max_merit = 123456789; // Initialize to an impossible value
     		halos_mpi[ihalo].id_max_merit = 123456789; 
	} else {     
		halos_mpi[ihalo].haloid_max_merit = halos[isimu1][ihalo].mtree[0].haloid[1]; 
     		halos_mpi[ihalo].id_max_merit = halos[isimu1][ihalo].mtree[0].id[1]; 
	}
  }

  return(0);
}

int MPI_Swaphalos_mpi(int isimu1)
{
	
	HALO_MPIptr halos_mpiRecvBuffer=NULL;
	size_t size_HALO_MPI = sizeof(HALO_MPI);
	
	uint64_t nHalosLocTmp=0;

	nHalosLocTmp = nHalos[isimu1];

	//  free(halos_mpiSendBuffer[isimu1]);
	 // halos_mpiSendBuffer[isimu1] = (HALO_MPIptr) calloc(nHalos[isimu1], sizeof(HALO_MPI));
	//  memcpy(&halos_mpiSendBuffer[isimu1]->haloid, &halos_mpi[isimu1]->haloid, nHalos[isimu1] * sizeof(HALO_MPI));

#ifdef DEBUG_MPI
	fprintf(stderr, "SwapMPI, Pre-send, nHalos=%"PRIu64", task=%d\n", nHalos[isimu1], LocTask);
#endif


 	  MPI_Sendrecv(&nHalos[isimu1], sizeof(uint64_t), MPI_BYTE, RecvTask, 0,
                &nHalosRecvBuffer[isimu1],  sizeof(uint64_t), MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

//	  MPI_Sendrecv(&nHalosLocTmp, sizeof(uint64_t), MPI_BYTE, RecvTask, 0,
 //            &nHalosRecvBuffer[isimu1], sizeof(uint64_t), MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

 // 	  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG_MPI
	fprintf(stderr, "SwapMPI - task %d sent %"PRIu64" halos to task=%d, received %"PRIu64" halos from task %d, TMP=%"PRIu64". \n ", 
		LocTask, nHalos[isimu1], RecvTask, nHalosBuffer[isimu1], SendTask, nHalosLocTmp);
#endif

	  halos_mpiRecvBuffer = (HALO_MPIptr) calloc(nHalosRecvBuffer[isimu1], size_HALO_MPI);

	  MPI_Sendrecv(halos_mpi, nHalos[isimu1] * size_HALO_MPI, MPI_BYTE, RecvTask, 0, 
	     halos_mpiRecvBuffer, nHalosRecvBuffer[isimu1] * size_HALO_MPI, MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

	  nHalos[isimu1] = nHalosRecvBuffer[isimu1];

#ifdef DEBUG_MPI
	fprintf(stderr, "SwapMPI, freeing up buffers... nHalos=%"PRIu64", task=%d\n", nHalos[isimu1], LocTask);
#endif

	  halos_mpi = (HALO_MPIptr) calloc(nHalos[isimu1], size_HALO_MPI); 	// TODO FIX ALL THE FREEs
	  memcpy(halos_mpi, halos_mpiRecvBuffer, nHalos[isimu1] * size_HALO_MPI);

	if(halos_mpiRecvBuffer != NULL) free(halos_mpiRecvBuffer);

	return(1);
}



int order_by_merit(int isimu0)
{
  MTREE    *mtree=NULL;
  uint64_t  ihalo, icroco, npart[2], jcroco, ncroco, common;
  double *merit=NULL;
  long unsigned *idx=NULL;

  npart[0] = 0;
  npart[1] = 0;

#ifdef DEBUG_MPI
	fprintf(stderr, "OrderByMerit on Task=%d.\n", LocTask);
#endif

    /* Order haloes in mtree by merit function */
    for(ihalo=0; ihalo<nHalos[isimu0]; ihalo++) 
    {
      ncroco = halos[isimu0][ihalo].global_ncroco;

      if(ncroco > 0)
      {

        mtree  = (MTREEptr)        calloc(ncroco, sizeof(MTREE));
        merit  = (double *) calloc(ncroco, sizeof(double));
        idx    = (long unsigned *) calloc(ncroco, sizeof(long unsigned));

	for(icroco=0; icroco<ncroco; icroco++)
        {
          common = halos[isimu0][ihalo].mtree[icroco].common;
          npart[0] = halos[isimu0][ihalo].mtree[icroco].npart[0]; 
          npart[1] = halos[isimu0][ihalo].mtree[icroco].npart[1]; 

          mtree[icroco].id[0]     = halos[isimu0][ihalo].mtree[icroco].id[0];
          mtree[icroco].haloid[0] = halos[isimu0][ihalo].mtree[icroco].haloid[0];
          mtree[icroco].npart[0]  = halos[isimu0][ihalo].mtree[icroco].npart[0];
          mtree[icroco].common    = halos[isimu0][ihalo].mtree[icroco].common;
          mtree[icroco].id[1]     = halos[isimu0][ihalo].mtree[icroco].id[1];
          mtree[icroco].haloid[1] = halos[isimu0][ihalo].mtree[icroco].haloid[1];
          mtree[icroco].npart[1]  = halos[isimu0][ihalo].mtree[icroco].npart[1];

          merit[icroco] = pow2((double)common)/((double)npart[0]*(double)npart[1]);
	} 

	    /* order by merit function */
  	 //   indexx((long unsigned)ncroco, merit-1, idx-1);	TODO fix this

	/* Free the mtree of halos and realloc it from scratch so we can store the tree in the right order (NO NEED TO DO IT!)*/
//	free(halos[isimu0][ihalo].mtree);
//	halos[isimu0][ihalo].mtree = (MTREEptr) calloc(ncroco, sizeof(MTREE));

	/* Now gather again the */
        for(jcroco=0; jcroco<ncroco; jcroco++)
        {
          icroco = idx[ncroco-1-jcroco]-1;

	  halos[isimu0][ihalo].mtree[jcroco].id[0] = mtree[icroco].id[0];
	  halos[isimu0][ihalo].mtree[jcroco].haloid[0] = mtree[icroco].haloid[0];
	  halos[isimu0][ihalo].mtree[jcroco].npart[0] = mtree[icroco].npart[0];
	  halos[isimu0][ihalo].mtree[jcroco].common = mtree[icroco].common;
	  halos[isimu0][ihalo].mtree[jcroco].id[1] = mtree[icroco].id[1];
	  halos[isimu0][ihalo].mtree[jcroco].haloid[1] = mtree[icroco].haloid[1];
	  halos[isimu0][ihalo].mtree[jcroco].npart[1] = mtree[icroco].npart[1];
        }
      }

    }
	
#ifdef DEBUG_MPI
	fprintf(stderr, "Free() OrderByMerit on Task=%d.\n", LocTask);
#endif
	//if(mtree != NULL)  free(mtree);	// TODO DON'T DO THESE FREEs HERE
	//if(idx != NULL)  free(idx);
 	//if(merit != NULL) free(merit);

  return(0);
}

#ifdef MTREE_BOTH_WAYS
/* Once the connection has been established backwards and forwards, determine the unique descendant of each halo */
int clean_connection(uint64_t ihalo, int isimu0, int isimu1, int iloop)
{
  uint64_t jhalo, khalo, icroco, loc_croco;
  
	// On the first loop this is just zero
	  loc_croco = ncroco_simu[ihalo];

     /* first make sure that the halo DOES have a connection */
     for(icroco=0; icroco<halos[isimu0][ihalo].global_ncroco; icroco++) 
     {
	/* this ID was local, however, we did not change the order of haloes in the halos_mpi files
         * so this should still identify the local halo position in the right chunk
         */
          jhalo = halos[isimu0][ihalo].mtree[icroco].id[1];
          khalo = halos[isimu0][ihalo].haloid;

	/* jhalo is the local id of the halo in the respective chunk. If jhalo is larger than nHalos[1] this means that this halo was found in a different 
	 * chunk so we skip the comparison altogether.
	 **/
	if(jhalo < nHalos[isimu1])
        if(khalo == halos_mpi[jhalo].haloid_max_merit) 
	{
			// TODO THIS IS JUST A MANUAL TWEAK - TO PREVENT HALO SPLITTING !
	      if((float) halos[isimu0][ihalo].mtree[icroco].common / (float) halos[isimu0][ihalo].mtree[icroco].npart[1] > 0.05 )
	      {
	    /* copy the new mtree accoring to the max merit */
	       ncroco_simu[ihalo]++;
	       loc_croco = ncroco_simu[ihalo];
	       mtree_tmp[ihalo] = (MTREE *) realloc(mtree_tmp[ihalo], (loc_croco + 1) * sizeof(MTREE));
	       mtree_tmp[ihalo][loc_croco-1].id[0]     = halos[isimu0][ihalo].mtree[icroco].id[0];
               mtree_tmp[ihalo][loc_croco-1].haloid[0] = halos[isimu0][ihalo].mtree[icroco].haloid[0];
	       mtree_tmp[ihalo][loc_croco-1].npart[0]  = halos[isimu0][ihalo].mtree[icroco].npart[0];
	       mtree_tmp[ihalo][loc_croco-1].common    = halos[isimu0][ihalo].mtree[icroco].common;
	       mtree_tmp[ihalo][loc_croco-1].id[1]     = halos[isimu0][ihalo].mtree[icroco].id[1];
	       mtree_tmp[ihalo][loc_croco-1].haloid[1] = halos[isimu0][ihalo].mtree[icroco].haloid[1];
	       mtree_tmp[ihalo][loc_croco-1].npart[1]  = halos[isimu0][ihalo].mtree[icroco].npart[1];
	       jhalo = 0;
	      }
	}
        //fprintf(stderr, "CleanConnection %"PRIu64", %"PRIu64", %"PRIu64", %"PRIu64" , %"PRIu64".\n", icroco, ihalo, iloop, nHalos[isimu0], nHalos[isimu1]);
     }

//        fprintf(stderr, "CleanConnection %"PRIu64" , %"PRIu64".\n", ihalo, iloop);
  /* replace halos[isimu0][ihalo].mtree[] with new structure array on the last loop */
  if (iloop == TotTask-1)
  {

    if(halos[isimu0][ihalo].mtree != NULL)
    {
       free(halos[isimu0][ihalo].mtree);
         halos[isimu0][ihalo].mtree = NULL;
    }

 //       fprintf(stderr, "CleanConnection %"PRIu64" , %"PRIu64".\n", ihalo, iloop);
    halos[isimu0][ihalo].global_ncroco = ncroco_simu[ihalo];
    halos[isimu0][ihalo].mtree = mtree_tmp[ihalo];

//fprintf(stderr, "(%d/%d) Clean connection, halo %llu ID = %llu has %llu/%llu cross correlations\n", 
//	iloop, TotTask, ihalo, halos[isimu0][ihalo].haloid, ncroco_simu[ihalo], halos[isimu0][ihalo].global_ncroco);
  }
  return 1;
}
#endif

/*==================================================================================================
 * max_merit
 *==================================================================================================*/
uint64_t max_merit(uint64_t jhalo)
{
  /* mtree[] is ordered by merit and hence we only need to check the first entry */
  if(halos_mpi[jhalo].ncroco > 0) {
    return(halos_mpi[jhalo].haloid_max_merit);
  }
  else {
#ifdef DEBUG_MPI
    fprintf(stderr,"jhalo=%ld in does not point to anywhere!?\n",jhalo);
#endif
    return(0);
  }
}

  /* Reallocate halos after having swapped particles among tasks */
int alloc_halos(int isimu)
{
  uint64_t ihalo=0;
#ifdef DEBUG_MPI
	fprintf(stderr, "Alloc %"PRIu64" halos for sim %d on task %d\n", nHalos[isimu], isimu, LocTask); 
#endif

/*
  //int n_halos_left = sizeof(halos[isimu])/sizeof(HALO);
//	fprintf(stderr, "Nhalos left %d \n", n_halos_left);

  if (halos[isimu] != NULL) 
  {
      for (ihalo=0; ihalo<nHalos[isimu]; ihalo++) 
      {
	fprintf(stderr, "Task %d, simu %d ihalo %d halo %"PRIu64" npart %"PRIu64"  ncroco %"PRIu64"\n",
		LocTask, isimu, ihalo, halos[isimu][ihalo].haloid, halos[isimu][ihalo].npart, halos[isimu][ihalo].ncroco);	
	if (halos[isimu][ihalo].npart > 0)
  	{
 	 //	if (halos[isimu][ihalo].Pid != NULL) //(uint64_t *) malloc(sizeof(uint64_t));
	//    free(halos[isimu][ihalo].Pid); 
	//    free(halos[isimu][ihalo].mtree); //(uint64_t *) malloc(sizeof(uint64_t));
	}
//	if (halos[isimu][ihalo].mtree != NULL) //(MTREEptr) malloc(sizeof(MTREE));
//	    free(halos[isimu][ihalo].mtree); //(MTREEptr) malloc(sizeof(MTREE));
 // 	if (halos[isimu][ihalo].Pid != NULL) //(uint64_t *) malloc(sizeof(uint64_t));
      }
  //    free(halos[isimu]);
  }
*/	// TODO this was a debug loop !!!

  halos[isimu] = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));

	/* Initialize pointers */ 
  for (ihalo=0; ihalo<nHalos[isimu]; ihalo++) 
  {
  	halos[isimu][ihalo].ncroco = 0;
  	halos[isimu][ihalo].npart = 0;
  	halos[isimu][ihalo].global_ncroco = 0;
  	halos[isimu][ihalo].mtree = NULL; 
  	halos[isimu][ihalo].Pid = NULL; 
  }

 return(1);
}


/* Realloc the haloes struct array for the newly read halos_tmp */
int add_halos(int ifile, int isimu)
{
  uint64_t nparts, ihalo, min_halo;
  size_t sizeHalo;
  
  if(ifile == 0) /* init the halo struct when reading the first file */ 
  {  
     min_halo = 0;
     nHalos[isimu] = nHalosTmp[isimu];
     totHaloSize = totHaloSizeTmp;
     alloc_halos(isimu);
  } else { /* realloc and make room for the new halos */
     min_halo = nHalos[isimu];
     nHalos[isimu] += nHalosTmp[isimu];
     totHaloSize += totHaloSizeTmp;
     halos[isimu] = (HALOptr) realloc(halos[isimu], ( (nHalos[isimu])) * sizeof(HALOS)); //TODO nHalos+1 ??
  }

  /* copy the halos into the old structure and free up memory */

  for(ihalo=min_halo; ihalo<nHalos[isimu]; ihalo++)
  {
     nparts = halos_tmp[isimu][ihalo-min_halo].npart;
     halos[isimu][ihalo].npart = halos_tmp[isimu][ihalo-min_halo].npart;
     halos[isimu][ihalo].haloid = halos_tmp[isimu][ihalo-min_halo].haloid;

	/* Initialize mtree and connections */
     halos[isimu][ihalo].mtree = NULL;
     halos[isimu][ihalo].global_ncroco = 0;

     memcpy(halos[isimu][ihalo].Xc, halos_tmp[isimu][ihalo-min_halo].Xc, 3 * sizeof(float));

     sizeHalo = nparts * sizeof(uint64_t);

     halos[isimu][ihalo].Pid = (uint64_t *) malloc(sizeHalo);

     memcpy(halos[isimu][ihalo].Pid, halos_tmp[isimu][ihalo-min_halo].Pid, sizeHalo);   

     /* free temp memory as the particles are copied */
     free(halos_tmp[isimu][ihalo-min_halo].Pid);
  }

  /* free the temporary structures */	
  free(halos_tmp[isimu]);

  nHalosTmp[isimu]=0;
  totHaloSizeTmp=0;

  return(1);
} 

  /* Clean up the halos */
int free_halos(int isimu)
{
  uint64_t ihalo;

  for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
  {
    if(halos[isimu][ihalo].Pid != NULL)
      free(halos[isimu][ihalo].Pid);

    if(halos[isimu][ihalo].mtree != NULL)
      free(halos[isimu][ihalo].mtree);
  }

 free(halos[isimu]);

 return(1);
}


/* custom linear time algorithm */
void intersection(int isimu0, int isimu1, uint64_t ihalo, uint64_t khalo, uint64_t *common)
{
	uint64_t ipart=0, kpart=0, apart=0, bpart=0;

	while (ipart < halos[isimu0][ihalo].npart && kpart < halos[isimu1][khalo].npart)
	{
		apart = halos[isimu0][ihalo].Pid[ipart];
		bpart = halos[isimu1][khalo].Pid[kpart];

		if (apart == bpart)
		{
			common[khalo]++;
			ipart++;
			kpart++;
		}
		else if (apart < bpart)
		{
			ipart++;
		}
		else
		{
			kpart++;
		}
	}
    
}

/* This function is called by qsort() and bsearch() needed for the matching of the particles */
//uint64_t cmpfunc(const void * a, const void * b)
uint64_t cmpfunc(uint64_t *a, uint64_t *b)
//int cmpfunc(const void * a, const void * b)
{
   return ( *(uint64_t*) a - *(uint64_t *) b  ) ;
//   return ( *(int*)a - *(int*)b );
}

uint64_t cmpdbl(double *a, double *b)
//int cmpfunc(const void * a, const void * b)
{
   return ( *(double *) a - *(double *) b  ) ;
//   return ( *(int*)a - *(int*)b );
}


int bubblesort(uint64_t *array, uint64_t size)
{
	uint64_t x=0, y=0, temp=0;

	for(uint64_t x=0; x<size; x++)
	{
		for(uint64_t y=0; y<size-1; y++)
		{

			if(array[y]>array[y+1])
			{

				temp = array[y+1];
				array[y+1] = array[y];
				array[y] = temp;

			}

		}

	}

  return 0;
}


int selectionsort(uint64_t *array, uint64_t size)
{
	uint64_t x=0, y=0, temp=0, index_of_min=0;

	for(x=0; x<size; x++)
	{

		index_of_min = x;

		for(y=x; y<size; y++)
		{

			if(array[index_of_min]>array[y])
			{

				index_of_min = y;

			}

		}

		temp = array[x];
		array[x] = array[index_of_min];
		array[index_of_min] = temp;
	}

	return 0;
}

/* construct a unique haloid based upon line number (=ihalo) and MPI_rank */
uint64_t constructHaloId(uint64_t ihalo)
{
  uint64_t haloid;
  
  // we allow for 64-48=16 bits for the MPI_rank (up to 65536 MPI tasks)
  haloid = (uint64_t) LocTask;
  haloid = haloid << 48;

  // simply overlay line number now
  haloid = haloid | ihalo;
  
  return(haloid);
}

/* Debug routines to check the content of the structures sent across tasks */
void check_halos(int isimu)
{
  unsigned int i=0, jmax=0, ntot;
  ntot = nHalos[isimu];
   jmax = 15;

    for(i=0; i<jmax; i++)
    {
      fprintf(stderr,"Local halo %d on task=%d has ID %"PRIu64" and %"PRIu64" parts Pid = %"PRIu64"\n", i, LocTask, halos[isimu][i].haloid, 
	halos[isimu][i].npart, halos[isimu][i].Pid[0]);
        //for(j=0; j<3; j++)
          // fprintf(stdout, "\t\tpart_id=%"PRIu64"\n", halos[isimu][i].Pid[j]);
    }

    for(i=ntot-jmax; i<ntot; i++)
    {
      fprintf(stderr,"Local halo %d on task=%d has ID %"PRIu64" and %"PRIu64" parts Pid = %"PRIu64"\n", i, LocTask, halos[isimu][i].haloid, 
	halos[isimu][i].npart, halos[isimu][i].Pid[0]);
       // for(j=0; j<3; j++)
         //  fprintf(stdout, "\t\tpart_id=%"PRIu64"\n", halos[isimu][i].Pid[j]);
    }

}


#ifdef DEBUG_LOG
void dump_log_halo(int isimu, int ihalo)
{
  uint64_t i=0, ntot, haloid;
  char log_name[200];
  FILE *log_halo;

  ntot = halos[isimu][ihalo].npart;
  haloid = halos[isimu][ihalo].haloid;

  sprintf(log_name, "/home/carlesi/MERGER_TREE/TEST/haloID_%03d.task_%02d.halo", ihalo, LocTask);
  log_halo = fopen(log_name, "w");

	//fprintf(stderr, "Printing halo %03d on task %d to %s\n", ihalo, LocTask, log_name);

	fprintf(log_halo, "# NPART=%"PRIu64"\t ID=%"PRIu64"\n", 
		halos[isimu][ihalo].npart,
		halos[isimu][ihalo].haloid);

    for(i=0; i<ntot; i++)
	fprintf(log_halo, "%"PRIu64"\n", halos[isimu][ihalo].Pid[i],
 
  fclose(log_halo);
}
#endif /* DEBUG_LOG */
