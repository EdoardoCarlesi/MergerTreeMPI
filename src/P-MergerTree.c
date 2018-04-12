 /*  MergerTree:   Merger Tree AHF_particles files
 *
 *
 *  input:    - how often to perform
 *            - 2x _particles files
 *
 *  output:   - 1x _mtree file
 *
 *
 * it is checked what halos in file2 make up the halos in file1, i.e.
 *
 *   file1   file2
 *
 *    0        0
 *    0       17
 *    0       31    -> halo #0 in file1 shares particles with halos #0,17,31 in file2
 *    1        2
 *    1       12
 *    1        4    -> halo #1 in file1 shares particles with halos #2,12,4  in file2
 *       etc.
 * TODO:
 * 	- Impose the HALO_MIN_PART already when reading particles in (save memory)
 *	- Use halos_mpi struct to send Xc[] - halo positions - and to read them in from AHF_halos
 *	- Check the periodic conditions algorithm (can be improved perhaps?)
 * 	- Each halo struct might just contain a pointer to halos_mpi instead of storing all the informations
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "param.h"
#include "tdef.h"
//#include "common.h"

#include "ptree_io.h"
#include "ptree_vars.h"
#include "ptree_funcs.h"

#include "libutility/utility.h"

#include <mpi.h>
//#include <omp.h>
/*==================================================================================================
 * main:
 *
 *       simply a wrapper for successive calls to create_mtree()
 *
 *==================================================================================================*/
int main(int argv, char **argc)
{
  int    ifile, jfile, jchunk, count, nFiles;   // nFiles is the number of snapshots (composed of several chunks)
  int 	 tmp_pos_mem=0, tmp_par_mem=0;
  char   *tempDir;		// Temporary files are stored here - to be cleaned at the end of the run
  char   *outList;		// List of snapshots file numbers, _??? format
  char   *partList;		// This file contains all the files to be submitted to task 0 
  char   *haloList;		// This file contains all the files to be submitted to task 0 
  char   *prefixOut;

  uint64_t ihalo=0;
  time_t   elapsed = (time_t)0, total = (time_t) 0;

  char   **outSuffix=NULL; 		// Suffix numbers to the particle files
  char   ***locPartFile=NULL;		// Each task stores _particle urls here
  char   ***locHaloFile=NULL;		// Each task stores _halo urls here, to read halo positions
  char   **locOutFile=NULL;		// Each task will dump to this file

  count = 1 ;
  
  /* set some global variables common to all MPI tasks */
  NReadTask = atoi(argc[count++]); 
  nFiles = atoi(argc[count++]);	  
  prefixOut = argc[count++];
  tempDir = argc[count++];
  outList = argc[count++];
  partList = argc[count++];	   
  haloList = argc[count++];	   
  BoxSize = atof(argc[count++]);	   
  Mparticle = atof(argc[count++]);

#ifdef DEBUG_INPUT
  /* check that the global variables have been read correctly */
  for(ifile=0; ifile<8; ifile++)
      fprintf(stderr, "argc[%d]=%s\n", ifile+1, argc[ifile+1]);
#endif

  /* general variables are set, now init MPI parallel part */
  MPI_Init(&argv, &argc);
  MPI_Comm_rank(MPI_COMM_WORLD, &LocTask); 
  MPI_Comm_size(MPI_COMM_WORLD, &TotTask); 

  /* RecvTask is recieving from LocTask and SendTask is sending to LocTask */
  SendTask = (LocTask+TotTask-1) % (TotTask); 
  RecvTask = (LocTask+1) % (TotTask);
  DeltaTask = TotTask - NReadTask; 

#ifdef USE_HALO_MIN_PART
  /* Send a warning message about ignoring smaller haloes */	
  if(LocTask == 0)
	fprintf(stderr, "Ignoring mergers for haloes with less than %d particles\n", HALO_MIN_PART);
#endif

  /* This is useful if you have more file-chunks than tasks */
  if(NReadTask > TotTask)
  {
     filesPerTask = (int) NReadTask / TotTask;  
     extraFilesPerTask = (int) NReadTask % TotTask;  

     if(LocTask < extraFilesPerTask)
       filesPerTask++; 
  }
  else
  {
     filesPerTask = 1;
  }

  /* allocate memory for locOutFiles, each of size MAXSTRING */
  locOutFile  = (char **) calloc(nFiles, sizeof(char *));
  outSuffix = (char **) calloc(nFiles, sizeof(char *));

  for(ifile=0; ifile<nFiles; ifile++)
    locOutFile[ifile]  = (char *) calloc(MAXSTRING, sizeof(char));

     /* each task has an otput file name */
     assign_output_files_names(outList, outSuffix, nFiles);

#ifdef DEBUG_MPI
  fprintf(stderr, "Task=%d is sending to and Task=%d recieving from LocTask=%d\n", SendTask, RecvTask, LocTask);
#endif

  if(LocTask == 0)
     fprintf(stderr, "MergerTree MPI mode, reading %d snapshots per task on %d tasks on a total of %d cpus.\n", 
		nFiles, NReadTask, TotTask);
  
  /* read the first file into memory */
  if(LocTask == 0)
    fprintf(stderr,"\nStartup:\n");

  if(LocTask < NReadTask)
  {
     /* alloc memory for local urls storage */
     locPartFile = (char ***) calloc(filesPerTask, sizeof(char **));
     locHaloFile = (char ***) calloc(filesPerTask, sizeof(char **));

     for(ifile=0; ifile<filesPerTask; ifile++)
     {
       locPartFile[ifile] = (char **) calloc(nFiles, sizeof(char *));
       locHaloFile[ifile] = (char **) calloc(nFiles, sizeof(char *));
     }

     /* now assign a list of urls to every task and alloc local memory */
     assign_input_files_to_tasks(partList, haloList, tempDir, locPartFile, locHaloFile, nFiles);
 
     /* read particles now reads into the tmp file which is reallocated in the main halos struct */
     for(ifile=0; ifile<filesPerTask; ifile++)
     {
        tmp_par_mem = read_particles(locPartFile[ifile][0], 0);
        tmp_pos_mem =read_positions(locHaloFile[ifile][0], 0);
        add_halos(ifile, 0);
     }
  }

    /*  If TotTask > NReadTask then redistribute particles across different tasks */
    if(DeltaTask > 0) 
      load_balance(0);  

    /* This is needed for faster shared particles comparison */
    order_halo_ids(0);

  if(LocTask == 0)
    fprintf(stderr,"\n");

  /* now loop over all file chunks to be read by each task */
  for(ifile=0; ifile<nFiles-1; ifile++)
  {
    sprintf(locOutFile[ifile], "%s%s%s.%04d", prefixOut, outSuffix[ifile], outSuffix[ifile+1], LocTask); 

    /* every task reads the next file into memory - the next file should be a chunk of _particle files
       at a different redshift */
  if(LocTask < NReadTask)
    for(jfile=0; jfile<filesPerTask; jfile++) /* this works if each task reads in more than one file */
    {
       /* be verbose */
       if(LocTask == 0 && jfile == 0)
         fprintf(stderr,"Correlating '%s' to '%s'\n           -> writing to '%s'\n",
           locPartFile[jfile][ifile],locPartFile[jfile][1+ifile],locOutFile[ifile]);
#ifdef DEBUG_MPI
  if(LocTask == 0)
       fprintf(stderr, "Loop=%d), jfile=%d, filename=%s\n", ifile, jfile, locPartFile[jfile][1+ifile]);
#endif
       read_particles(locPartFile[jfile][1+ifile], 1);
       read_positions(locHaloFile[jfile][1+ifile], 1);
       add_halos(jfile, 1);
    }

    /* if NTask > N files to be read then scatter the data through the tasks */
    if(DeltaTask > 0)   
       load_balance(1);  

    /* This is needed for faster shared particles comparison */
    order_halo_ids(1);

    /* now swap the files across all the tasks and look for correlations */
    for(jchunk=0; jchunk<TotTask; jchunk++)
    {
      elapsed -= time(NULL);

      /* cross correlate locHaloFile[i] to locHaloFile[i+1] */
      cross_correlation(jchunk, 0, 1);

      elapsed += time(NULL);
      total += elapsed;

      if(LocTask == 0 && jchunk == TotTask-1)
        fprintf(stderr, "\nCross correlation completed step %d/%d of file %d in %ld sec, total elapsed time is %ld s.\n", 
	        jchunk+1, TotTask, jfile, elapsed, total);

      MPI_Barrier(MPI_COMM_WORLD);

      /* Now that the comparison has been done, swap the halos structs through tasks. No need to swap on the last step, the loop is over */
      MPI_Swaphalos(1);

      elapsed = (time_t) 0;

#ifdef DEBUG_MPI
      fprintf(stderr, "Task= %d has received from task = %d a total of %"PRIu64" haloes\n",
		LocTask, SendTask, nHalos[1]);
      fprintf(stderr, "Total elapsed time on task=%d is %ld sec.\n", LocTask, total);
#endif

    } /* End swapping data */

#ifndef MTREE_BOTH_WAYS

	write_mtree(0, locOutFile[ifile]);  /* Write the merger graph - not the tree */

#else
      if(LocTask == 0)
	fprintf(stderr, "Starting forward cross correlation.\n");

    /* now swap the files across all the tasks and look for correlations 
     * this time isimu=1 is fixed and isimu=0 is swapped through the tasks
     * however, we keep the original halo[0] struct intact on each task, to "save" its merger tree
     */

    /* copy the particles and halos into new structures */
    nHalos[2] = nHalos[0];
	
    alloc_halos(2);

    /* First we copy the 0 halo structs into 2, since we will then swap these ones */
    copy_halos(0, 2);

    order_halo_ids(2);

    for(jchunk=0; jchunk<TotTask; jchunk++)
    {
      elapsed -= time(NULL);

      /* cross correlate locHaloFile[i+1] to locHaloFile[i] */
      cross_correlation(jchunk, 1, 2);

      elapsed += time(NULL);
      total += elapsed;

#ifdef DEBUG_MPI
      if(LocTask == 0)
        fprintf(stderr, "Forward cross correlation completed step %d/%d of file %d in %4.2f sec,total elapsed time is %4.2f s.\n", 
		jchunk+1, TotTask, jfile, (double)elapsed/CLOCKS_PER_SEC, (double)total);
#endif

      elapsed = (time_t) 0;

      if(LocTask == 0 && jchunk == 0)
	fprintf(stderr, "Swapping haloes across tasks...");

	/* swap the halos across tasks - don't do it on the last step */
      //if(jchunk != TotTask-1) 
	MPI_Swaphalos(2);

      if(LocTask == 0 && jchunk == TotTask-1)
	fprintf(stderr, " done.\n");

#ifdef DEBUG_MPI
	fprintf(stderr, "ForwardCorrelation: Task=%d has received %"PRIu64" halos from task=%d\n",
	 		LocTask, nHalosRecvBuffer[2], SendTask);
#endif
    } /* End swapping data for forward cross_correlation */

	/* we don't need the "temporary" halos[2] anymore */
	if(halos[2] != NULL) free_halos(2);

	/* now order by merit function the forward correlations for the mtree of the local 
	 * halo[1] of simu1 and store them into a temporary array of structures (since when
	 * we will swap the pointer insite the HALO struct we won't pass the whole mtree)
	 */
	order_by_merit(1);
	jchunk = 0;

	halos_mpi = (HALO_MPIptr) calloc(nHalos[1], sizeof(HALO_MPI));

      if(LocTask == 0)
	fprintf(stderr, "Copy haloes into MPI structures...");
	copy_halos_to_halos_mpi(1);

      if(LocTask == 0)
	fprintf(stderr, " done.\n");

	/* this variable stores the cross_connections that uniquely connect a progenitor to a descendant halo ON EACH TASK */
	ncroco_simu = (uint64_t *) calloc(nHalos[0], sizeof(uint64_t));

	/* Initialize the temp mtree struct to store all the data while clearing the connection, then it will be copied to the halos[][].mtree */
     	mtree_tmp = (MTREEptr *) calloc(nHalos[0], sizeof(MTREEptr));

	for (ihalo=0; ihalo<nHalos[0]; ihalo++)
	{
     		mtree_tmp[ihalo] = NULL; //(MTREEptr) calloc(1, sizeof(MTREE));
     		ncroco_simu[ihalo] = 0;
	}

  /* now clean connections simu0->simu1 
   * we need to swap all halos_mpi in simu1 again for proper comparison...
   */
    for(jchunk=0; jchunk<TotTask; jchunk++)
    {
	/* First clean connection simu0->simu1 using local data */ 
#ifdef DEBUG_MPI
      if(LocTask == 0)
	fprintf(stderr, "clean_connection() 0->1 (%d/%d) on task=%d\n", jchunk+1, TotTask, LocTask);
#endif
 	  for(ihalo=0; ihalo<nHalos[0]; ihalo++)	
 	      clean_connection(ihalo, 0, 1, jchunk);

	  MPI_Swaphalos_mpi(1);

#ifdef DEBUG_MPI
	fprintf(stderr, "POSTSWAPMPI - %d) LocTask=%d, nhalos1=%"PRIu64", nhalos2=%"PRIu64" \n", jchunk, LocTask, nHalos[0], nHalos[1]);
      if(LocTask == 0)
	fprintf(stderr, " done.");
#endif

	/* if(mtree_tmp != NULL) 
	{
		for (ihalo=0; ihalo<nHalos[0]; ihalo++)
     		if (mtree_tmp[ihalo] != NULL) free(mtree_tmp[ihalo]); //(MTREEptr) calloc(1, sizeof(MTREE));
		free(mtree_tmp);
	}*/
  } /* End loop on jchunk of halos[1] */

     write_mtree(0, locOutFile[ifile]);  
	//free_halos(1);
	// CHECK ALL THESE FREEEs TODO
 // for (ihalo=0; ihalo<nHalos[1]; ihalo++) 
 // {    
	//if(mtree_tmp[ihalo] != NULL)
		//free(mtree_tmp[ihalo]);
		//mtree_tmp[ihalo] = NULL;
 // }
 //    if(mtree_tmp != NULL) free(mtree_tmp);
     free(ncroco_simu);		
     free(halos_mpi);

#endif /* MTREE_BOTH_WAYS*/

    /* be verbose */
    if(LocTask == 0)
      fprintf(stderr,"  o making file 1 the new file 0 ...");

    /* remove halo[0] structs from memory */
    free_halos(0);

    
    /* make HaloFile[i+1] the new HaloFile[i] */
    nHalos[0] = nHalos[1];
    alloc_halos(0);
    copy_halos(1,0);
    free_halos(1);

//    halos[0] = halos[1];

    /* be verbose */
    if(LocTask == 0)
    	fprintf(stderr," done\n");

  } /* for(nFiles) */

  /*==================================================================*
   *                             CLEANUP                              *
   *==================================================================*/
    if(LocTask == 0)
      fprintf(stderr,"\nCleaning up ...\n");
	
exit(0);
//  if (tmp_par_mem > 0 && tmp_pos_mem > 0)
  {
  /* free input filenames - set it as an extra option */
  if(LocTask < NReadTask)
  {
     for(jfile=0; jfile<filesPerTask; jfile++)
     {
        for(ifile=0; ifile<nFiles; ifile++)
	{
          if(locPartFile[jfile][ifile]) free(locPartFile[jfile][ifile]);
          if(locHaloFile[jfile][ifile]) free(locHaloFile[jfile][ifile]);
	}
        if(locPartFile[jfile]) free(locPartFile[jfile]);
        if(locHaloFile[jfile]) free(locHaloFile[jfile]);
     }
     if(locPartFile) free(locPartFile);
     if(locHaloFile) free(locHaloFile);
  }

  /* free output filenames */
  for(ifile=0; ifile<nFiles; ifile++)
    if(locOutFile[ifile])  free(locOutFile[ifile]);

  if(locOutFile)  free(locOutFile);

  MPI_Finalize();
  }
	

  if(LocTask == 0)
    printf("finished\n");

  return(0);
}


