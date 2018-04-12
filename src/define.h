#ifndef DEFINE_INCLUDED
#define DEFINE_INCLUDED

				/* Comparison Options */
#define MINCOMMON       15        		// we only cross-correlate haloes if they at least share MINCOMMON particles
//HALO_MAX_DISTHALO_MIN_PART		// If the option is switched on, we neglect all cross correlations for haloes with less particles. Haloes with more particles are correlated with haloes smaller that this number, so that minor mergers are included (but the final mtree catalogs will only feature larger haloes)
#define	USE_HALO_MIN_PART
#define HALO_MIN_PART	100		// This gets rid of lots of haloes, keeps mayor mergers only. However it keeps track of smaller merger to larger haloes


#define MAX_HALO_DIST 3000 	// Maximum c.o.m. distance between halos - if larger than this, we skip the comparison (kpc)
//#define DYN_MAX_DIST		// This option computes the radius of comparison of haloes dynamically, the base value is always the MAX_HALO_DIST, the base value is always the MAX_HALO_DIST - TODO check it again, it seems it is slowing down things...
//#define ONLY_USE_PTYPE 1
#define MTREE_BOTH_WAYS		// make sure that every halo has only one descendant


				/* Code Options*/
//#define IMPROVE_LB			// Trying to improve the loadbalance between tasks. Not extensively tested.
#ifdef IMPROVE_LB
#define LB_PART_FAC 1.0
#define LB_HALO_FAC 1.0
#endif
				/* Output Options */
#define NO_HEADER		// DO NOT PRINT headers for the files
//#define SUSSING2013                   // write _mtree in format used for Sussing Merger Trees 2013

				/* Extra Options */
//#define DEBUG_MPI
//#define DEBUG_LOG


#define VERSION 1.0
#define BUILD   91

#ifdef AHF2
  #undef  VERSION
  #define VERSION 2.0
  #undef  BUILD
  #define BUILD   0
#endif

#define TERMINATE_AMIGA {"terminateAMIGA"}
#define AHF_MINPART_GAS     10     /* min. number of gas for spin and shape calculation                                          */
#define AHF_MINPART_STARS   10     /* min. number of stars for spin and shape calculation                                        */
#define AHF_MINPART_SHELL   10     /* min. number of particles in a profile shell for using AHFshellshape                        */
#define AHF_NBIN_MULTIPLIER 1      /* (integer value) increases the number of radial bins by this factor from the standard value */
#define AHF_HIRES_DM_WEIGHT 1.0
#define AHF_HOSTHALOLEVEL   1      /* first level to be considered as credible to spawn subbaloes                                */
#define AHF_HOSTSUBOVERLAP  0.5    /* how far should the subhalo have entered into the host                                      */
#define AHF_MIN_REF_OFFSET  0      /* offset for first refinement to be used by AHF                                              */
#define AHF_RISE            1.00   /* Rho > AHF_RISE*Rho_prev -> rising density                                                  */
#define AHF_SLOPE           0.99   /* outer halo profile at least like r^-AHF_SLOPE                                              */
#define AHF_MAXNRISE        2      /* try to catch variations in density                                                         */
#define AHF_Rmax_r2_NIGNORE 5      /* how many central particle to ignore with AHFparticle_RMax_r2 feature                       */
#define PGAS                0.0   /* identifier for gas particles; has to be exactly 0.0!!!!                                     */
#define PDM                -1.0   /* identifier for dm particles; whatever negative value                                        */
#define PSTAR              -4.0   /* identifier for star particles; whatever negative value                                      */
#define PDMbndry           -5.0   /* identifier for bndry particles; whatever negative value                                     */

#define MIN_NNODES	15


#endif

