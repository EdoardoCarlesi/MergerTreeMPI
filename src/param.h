#ifndef PARAM_INCLUDED
#define PARAM_INCLUDED

#include <float.h>

/*================================================================================
 * switch on/off various features of AMIGA not passed via DEFINEFLAGS in Makefile
 *================================================================================*/
#include "define.h"

/* define */

#define COUNT_VAR	50000	// Print the output every COUNT_VAR cycles

/*=============================================================================
 * some physical constants...
 *=============================================================================*/
#define Gyr       3.1558e16         /* [sec]                   */
#define Mpc       3.08567782e19     /* [km]                    */
#define H0        100.              /* [h*km]/[sec*Mpc]        */
#define rhoc0     2.7755397e11      /* [h^2*Msun]/[Mpc^3]      */
#define Grav      4.3006485e-9      /* [Mpc*km^2]/[Msun*sec^2] */
#define cH0	      2998.0		        /* c/H0 (in h^-1 Mpc)      */
#define kB_per_mp 0.825481286614E-2 /* [(km/sec)^2/K]          */
#define kBoltzman 6.9416792         /* [(km/sec)^2 Msun/K]     */
#define Msun      1.9891e30         /* [kg]                    */

/*============================================================================
 * handling output file names etc.
 *============================================================================*/
#define AMIGAHEADER  2048   /* maximum size (in bytes) for output file header   */
#define HEADERSTRING 256    /* no. of characters for header string in outfiles  */
#define HEADERSIZE  (HEADERSTRING*sizeof(char)+2*sizeof(long)+6*sizeof(int)+46*sizeof(double))
#define FILLHEADER  (AMIGAHEADER-HEADERSIZE)

#define bytes2GB  9.313225746154785e-10

#define MAXSTRING 1024	
#define MAXTIME 100

/** Number of bits used per dimension for the calculation of the Hilbert key */
#define BITS_PER_DIMENSION 21

/** The maximum amount of particles that can be send in on flush */
#define MAX_SEND_PARTICLES 1000000

/** Defines the verbosity level for the io_logging_* functions if
 * NEWSTARTUN is used. Depending on this value some messages might not
 * appear and hence this can be used to reduce the chatter produced by
 * the io_logging_* function. During the starup this will be passed to
 * the logging module (in startrun.c). The lower the number, the less
 * output will be produced. */
#define VERBOSITY 6

/*============================================================================
 * time stepping
 *============================================================================*/
#define NSTEPS       1000    /* number of (initial) steps for time stepping     */
#define CA_CRIT      0.15    /* restricts timestep due to da/a criterion    0.05*/ 
#define CELLFRAC_MAX 0.2     /* how far are particles allowed to move       0.2*/
#define CELLFRAC_MIN 0.05    /* how far should particles move at least      0.05*/
#define CF_MEAN      ((CELLFRAC_MAX+CELLFRAC_MIN)/2)

/*=============================================================================
 * finally some convenient abreviations...
 *=============================================================================*/
#define NDIM      3          /* DO NOT EVER TOUCH THIS NUMBER!  */
#define CRITMULTI 8.0
#define NP_RATIO  7.5
#ifdef DOUBLE 
#define ZERO         (1E-12)
#define MACHINE_ZERO (5E-16)
#else
#define ZERO         (1e-6)
#define MACHINE_ZERO (5e-16)
#endif
#define MZERO        (1e-10)  /* used when reading GADGET files */


/*=============================================================================
 * GADGET related parameters (only relevant for libio_serial.a !!!)
 *=============================================================================*/
#define GADGET_MUNIT 1.0e10    /* GADGET mass unit in Msol/h                   */
#ifdef GADGET_LUNIT_KPC
#define GADGET_LUNIT 1.0e-3
#else                          /* GADGET length unit in Mpc/h                  */
#define GADGET_LUNIT 1.0
#endif

/*============================================================================= 
 * no code is complete without defining these numbers ;-)
 *=============================================================================*/
#define PI    3.14159265358979323846264338
#define TWOPI 6.28318530717958647692528677
#define SQRT2 1.41421356237309504880168872

/*=============================================================================
 * The numbers 0, 1 and 2 are used to represent the X, Y and Z coordinates of
 * a three dimensional object (eg vector). To aid readability these numbers are
 * replaced by the names X, Y and Z when used individually.
 *=============================================================================*/
#define X 0     /* x-coord symbol */
#define Y 1     /* y-coord symbol */
#define Z 2     /* z-coord symbol */

/*=============================================================================
 * Boolean parameters
 *=============================================================================*/
#define YES   1
#define NO    0

#define ON    1
#define OFF   0

#define TRUE  1
#define FALSE 0

#endif

