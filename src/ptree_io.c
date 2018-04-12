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

/*==================================================================================================
 * assign_input_files_to_tasks:
 *
 * assigns to the local task the list of files which will have to be read.
 * For the moment, each task reads in one file per redshift.
 *
 *      partList  = a list with all the particle files which will be read by task 0
 *      tempDir	  = a temporary folder where all the lists of files corresponding to each task are stored
 *      NReadTask = number of processors actually reading (must be smaller or equal to the number of tasks) 
 *
 *==================================================================================================*/
int assign_input_files_to_tasks(char *partList, char *haloList, char *tempDir, char ***locPartFile, char ***locHaloFile, int nFiles)
{
  int ifile, jfile;
  FILE *locPartListFile=NULL;
  FILE *locHaloListFile=NULL;
  char command[MAXSTRING];
  char locPartList[MAXSTRING];
  char locHaloList[MAXSTRING];
  char dummy[MAXSTRING-1];
  char task[6];
  
  if(LocTask == 0)
    fprintf(stderr, "\n  o assigning files to each task from %s and %s\n", partList, haloList);

  //locPartList[MAXSTRING] = (char *) calloc(filesPerTask, sizeof(char));

  if(LocTask < NReadTask)
  {
    
    for(ifile=0; ifile<filesPerTask; ifile++)
    {
      sprintf(task, "%04d", LocTask + ifile * TotTask); 
      sprintf(locPartList, "%s/part.%s.tmp", tempDir, task);
      sprintf(command, "sed s/0000/%s/ <%s >%s", task, partList, locPartList);
      system(command);
      sprintf(locHaloList, "%s/halo.%s.tmp", tempDir, task);
      sprintf(command, "sed s/0000/%s/ <%s >%s", task, haloList, locHaloList);
      system(command);

      locPartListFile = fopen(locPartList, "r");

      for(jfile=0; jfile<nFiles; jfile++)
      {
        fgets(dummy, MAXSTRING-1, locPartListFile);
        locPartFile[ifile][jfile] = (char*) calloc(strlen(dummy)+5, sizeof(char));
        strcpy(locPartFile[ifile][jfile], dummy);
        locPartFile[ifile][jfile][strlen(dummy)-1]='\0';
#ifdef DEBUG_MPI
        fprintf(stderr, "Task=%d will read from file %s, list=%s\n", LocTask, locPartFile[ifile][jfile], locPartList);
#endif
      }
	/* Assing the halo files to each task */
      locHaloListFile = fopen(locHaloList, "r");

      for(jfile=0; jfile<nFiles; jfile++)
      {  
        fgets(dummy, MAXSTRING-1, locHaloListFile);
        locHaloFile[ifile][jfile] = (char*) calloc(strlen(dummy)+5, sizeof(char));
        strcpy(locHaloFile[ifile][jfile], dummy);
        locHaloFile[ifile][jfile][strlen(dummy)-1]='\0';
#ifdef DEBUG_MPI
        fprintf(stderr, "HALO Task=%d will read from file %s, list=%s\n", LocTask, locHaloFile[ifile][jfile], locHaloList);
#endif
      }

      /* close the temp file generated with the list of the files to be read */
      fclose(locPartListFile);
      locPartListFile = NULL;
    } /* ifile, loop on the file chunks that each task is loading in */
  } /* else task remains idle */

  if(LocTask == 0)
    fprintf(stderr, " done.\n");

  return(1);
}


/* 
 * This functions assign to each task its corresponding output file.
 */
int assign_output_files_names(char *outList, char **outSuffix, int nFiles)
{
  FILE *locOutListFile;
  int ifile=0;
  char dummy[MAXSTRING-1];

    locOutListFile = fopen(outList, "r");

      for(ifile=0; ifile<nFiles; ifile++)
      {
        fgets(dummy, MAXSTRING-1, locOutListFile);
        outSuffix[ifile] = (char*) calloc(strlen(dummy)+5, sizeof(char));
        strcpy(outSuffix[ifile], dummy);
        outSuffix[ifile][strlen(dummy)-1]='\0';
#ifdef DEBUG_MPI
        fprintf(stderr, "Task=%d will out to file %s\n", LocTask, outSuffix[ifile]);
#endif
      }

      /* close the temp file generated with the list of the files to be read */
      fclose(locOutListFile);

  return(1);
}


/*==================================================================================================
 * read_particles:
 *
 * read the file storing the particle IDs for each halo
 *
 *      nHalos = number of halos found in file
 *      Pid    = id's of all those particles
 *
 *==================================================================================================*/
int read_particles(char filename[MAXSTRING], int isimu)
{
  FILE     *fpin;
  char      line[MAXSTRING];
  int64_t   ihalo;
  uint64_t  nPartInHalo, nPartInUse, ipart, Pid, Ptype, haloid;
  time_t    elapsed = (time_t)0;

  elapsed -= time(NULL);

#ifdef DEBUG_MPI
    fprintf(stderr,"  o reading _particles file %s on task %d ...\n",filename, LocTask);
#else
  if(LocTask == 0)
    fprintf(stderr,"  o reading _particles file %s on task %d ...",filename, LocTask);
#endif

  fpin = fopen(filename,"r");
  if(fpin == NULL)
   {
    fprintf(stderr,"could not open file %s on task %d \nexiting!\n",filename, LocTask);
    exit(0);
   }
  
  /* reset all variables */
  nHalosTmp[isimu] = 0;
  ihalo         = -1; // -1 so it was before
  halos_tmp[isimu]  = NULL;

  /* this array keeps track of the size of each halo's dynamically allocated particle informations, Pid */
  totHaloSizeTmp = 0;

  /* get the first line from file */
  fgets(line,MAXSTRING,fpin);
 
  /* for AHF_particles files the first line is numGoodHalos which we can happily ignore */
  if(strncmp(line,"#",1) != 0 && sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)

    fgets(line,MAXSTRING,fpin);  

  do {
    if(strncmp(line,"#",1) != 0)
     {
      /* has a haloid been written */
      if(sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)
       {
        /* if not, just get the number of particles */
        sscanf(line,"%"SCNu64, &nPartInHalo);
        
        /* construct a meaningful haloid using linenumber and MPI_rank */
        haloid = constructHaloId((uint64_t)(ihalo+1)); // +1, because ihalo has not been incremented yet!

       }

      /* found yet another halo */
      ihalo++;
      nHalosTmp[isimu] += 1;
      halos_tmp[isimu]   = (HALOptr) realloc(halos_tmp[isimu], nHalosTmp[isimu]*sizeof(HALOS));
      
      /* halos_tmp[][].Pid will be incrementally filled using realloc() */
      halos_tmp[isimu][ihalo].Pid = NULL;

      nPartInUse = 0;

      for(ipart=0; ipart<nPartInHalo; ipart++)
       {
        /* read line containing ID and possibly some more information */
        fgets(line,MAXSTRING,fpin);
        
        /* check whether we are able to read a meaningful particle type, too */
        if(sscanf(line,"%"SCNu64" %"SCNu64, &Pid, &Ptype) == 1) {
          /* if not, set Ptype to 1 as this is the type we will use below */
          Ptype = 1;
         }
        
        // here we can restrict the cross-correlation to a ceratain sub-set of all particles
#ifdef ONLY_USE_PTYPE
        if(Ptype == ONLY_USE_PTYPE)
#endif
         {
          halos_tmp[isimu][ihalo].Pid    = (uint64_t *) realloc(halos_tmp[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(uint64_t));
          if(halos_tmp[isimu][ihalo].Pid == NULL) {
            fprintf(stderr,"read_particles: could not realloc() halos_tmp[%d][%ld].Pid for %"PRIu64"particles\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1));
            exit(-1);
          }

          halos_tmp[isimu][ihalo].Pid[nPartInUse] = Pid;
          nPartInUse++;
         }
       }
      
      /* store number of particles in halo */
      halos_tmp[isimu][ihalo].npart = nPartInUse;  

      /* add to the total number of particles on task, plus the haloid and haloindex */    
      totHaloSizeTmp += 2 * (nPartInUse + 1) * sizeof(uint64_t) + 3 * sizeof(float);
     }

  } while( fgets(line,MAXSTRING,fpin) != NULL);
 
  fclose(fpin);
  
  elapsed += time(NULL);
	// TODO THIS IS ADDED JUST BECAUSE
//	nHalosTmp[isimu]--;
  if(LocTask == 0)
   fprintf(stderr,"\n done in %4.2f sec.\n Added %"PRIu64" haloes, total temporary buffer is %zd MB.\n", 
		(double)elapsed, nHalosTmp[isimu], totHaloSize/1024/1024); 

  return (totHaloSize/1024/1024);
}

/* Read halo positions from the AHF_halos files. We only need halo positions for later comparison of halo distances, 
   to avoid unnecessary comparison between distant haloes */
int read_positions(char filename[MAXSTRING], int isimu)
{
  FILE     *fpin;
  char      line[MAXSTRING];
  int64_t   ihalo;
  uint64_t  ID_dummy, Host_dummy;
  int       Nsub_dummy, Npart_dummy;
  double    Mvir_dummy, Xc, Yc, Zc;
  time_t    elapsed = (time_t)0;

  elapsed -= time(NULL);

#ifdef DEBUG_MPI
    fprintf(stderr,"  o reading _halos file %s on task %d ...\n",filename, LocTask);
#else
  if(LocTask == 0)
    fprintf(stderr,"  o reading _halos file %s on task %d ...",filename, LocTask);
#endif

  fpin = fopen(filename,"r");
  if(fpin == NULL)
   {
    fprintf(stderr,"could not open file %s on task %d \nexiting!\n",filename, LocTask);
    exit(0);
   }

  ihalo = 0;

  /* get the first line from file */
  fgets(line,MAXSTRING,fpin);
 
  /* for AHF_halos files the first line is the header which we can happily ignore */
  do {
     if(strncmp(line,"#",1) != 0)
     {
       sscanf(line,"%"SCNu64" %"SCNu64" %d %lf %d %lf %lf %lf", 
	&ID_dummy, &Host_dummy, &Nsub_dummy, &Mvir_dummy, &Npart_dummy, &Xc, &Yc, &Zc);
       
          halos_tmp[isimu][ihalo].haloid = ID_dummy;
          halos_tmp[isimu][ihalo].Xc[0] = Xc;
          halos_tmp[isimu][ihalo].Xc[1] = Yc;
          halos_tmp[isimu][ihalo].Xc[2] = Zc;

	ihalo++;
    }
  } while( fgets(line,MAXSTRING,fpin) != NULL);
 
  fclose(fpin);
  
  elapsed += time(NULL);

  if(LocTask == 0)
   fprintf(stderr," done in %4.2f sec.\n On Task %d, buffer contains %"PRIu64" haloes, total size=%zd MB\n", 
		(double)elapsed, LocTask, nHalos[isimu], totHaloSize/1024/1024); 

  return (totHaloSize/1024/1024);
}


/*==================================================================================================
 * write_mtree
 *==================================================================================================*/
int write_mtree(int isimu0, char OutFile[MAXSTRING])
{
  uint64_t  ihalo, jhalo, icroco;
  FILE *fpout, *fpout_idx;
  char outname[MAXSTRING], outname_idx[MAXSTRING];
  time_t   elapsed = (time_t)0;

  order_by_merit(isimu0);

  elapsed -= time(NULL);
  
  if (LocTask == 0)
  fprintf(stderr,"  o writing cross-correlation for %"PRIu64" haloes on task %d...\n",nHalos[0], LocTask);
  
  sprintf(outname,"%s_mtree",OutFile);
  strcpy(outname_idx, outname);
  strcat(outname_idx, "_idx");
  
  fpout = fopen(outname,"w");
  if(fpout == NULL) {
    fprintf(stderr,"could not open file %s on task %d,\nexiting\n",outname, LocTask);
    exit(0);
  }
  
  fpout_idx = fopen(outname_idx,"w");
  if(fpout_idx == NULL) {
    fprintf(stderr,"could not open file %s on task %d,\nexiting\n",outname_idx, LocTask);
    exit(0);
  }
  
#ifdef SUSSING2013
  fprintf(fpout,"%"PRIu64"\n",nHalos[0]);
#else // SUSSING2013

#ifndef NO_HEADER
  fprintf(fpout,"#   HaloID(1)   HaloPart(2)  NumProgenitors(3)\n");
  fprintf(fpout,"#      SharedPart(1)    HaloID(2)   HaloPart(3)\n");
  fprintf(fpout_idx,"# HaloID(1) HaloID(2)\n");
#endif // NO_HEADER

#endif // SUSSING2013
  fflush(fpout);
  fflush(fpout_idx);

  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {

    if(DeltaTask > 0) /* in this case you had to load balance and the largest halos will be in the last positions */
      jhalo = nHalos[0] - ihalo -1;
    else
      jhalo = ihalo;

    if(halos[isimu0][jhalo].global_ncroco > 0) 
    //if(halos[isimu0][jhalo].ncroco > 0 && halos[isimu0][jhalo].global_ncroco > 0)	
    {
     fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", 
	halos[isimu0][jhalo].mtree[0].haloid[0], halos[isimu0][jhalo].mtree[0].haloid[1]);
      fflush(fpout_idx);
      
#ifdef SUSSING2013
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[isimu0][jhalo].haloid,
              halos[isimu0][jhalo].global_ncroco);
#else // SUSSING2013
      fprintf(fpout,"%"PRIu64"  %"PRIu64"  %"PRIu64"\n",
              halos[isimu0][jhalo].haloid,
              halos[isimu0][jhalo].npart,
              halos[isimu0][jhalo].global_ncroco);
#endif // SUSSING2013
      fflush(fpout);
      
      for(icroco=0; icroco<halos[isimu0][jhalo].global_ncroco; icroco++) {
#ifdef SUSSING2013
        fprintf(fpout,"%"PRIu64"\n",
                halos[isimu0][jhalo].mtree[icroco].haloid[1]);
#else // SUSSING2013
        fprintf(fpout,"  %"PRIu64"  %"PRIu64"  %"PRIu64"\n",
                halos[isimu0][jhalo].mtree[icroco].common,
                halos[isimu0][jhalo].mtree[icroco].haloid[1],
                halos[isimu0][jhalo].mtree[icroco].npart[1]);
#endif // SUSSING2013
        fflush(fpout);
      }
    }
#ifdef SUSSING2013
    else {
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[isimu0][jhalo].haloid,
              halos[isimu0][jhalo].ncroco);      
    }
#endif // SUSSING2013
  }
  
  /* close files */
  fclose(fpout);
  fclose(fpout_idx);

  elapsed += time(NULL);

   if(LocTask == 0)
    fprintf(stderr,"\n done in %4.2f sec.\n", (double)elapsed);
  return(1);
}
