/*  Last edited: Aug 28 15:08 1995 (esr) */

/* $Id: dotter.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */
/* 
   dotter.h - Public interface for dotter

   Only 3 parameters are mandatory, the rest can be set to NULL.
   A minimal call would look like:

   dotter(type, 0, 0, qseq, 0, 0, sseq, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

   NOTE: qseq and sseq must be messalloc'ed in the calling routine.  
   They are messfree'd by Dotter.
*/

#include "blxview.h"   /* For MSP struct */

Graph 
dotter (

	char  type,        /* Mandatory, one of { P, N, X } 
			      P -> Protein-Protein
			      N -> DNA-DNA
			      X -> DNA-Protein */
	
	char *opts,        /* Optional, may be NULL 
			      Various options for display features */

	char *queryname,   /* Optional, may be NULL 
			      Name of Horizontal sequence */
	
	char *queryseq,	   /* Mandatory, NULL terminated string
			      Horisontal sequence - messfree'd by Dotter */

	int   qoff,	   /* Optional, may be NULL
			      Coordinate offset of horisontal sequence */

	char *subjectname, /* Optional, may be NULL 
			      Name of vertical sequence */

	char *subjectseq,  /* Mandatory, NULL terminated string
			      vertical sequence - messfree'd by Dotter */

	int   soff,	   /* Optional, may be NULL 
			      Coordinate offset of horisontal sequence */

	int   qcenter,	   /* Optional, may be NULL 
			      Coordinate to centre horisontal sequence on */

	int   scenter,	   /* Optional, may be NULL 
			      Coordinate to centre horisontal sequence on */

	char *savefile,	   /* Optional, may be NULL 
			      Filename to save dotplot to. Invokes batch mode */

	char *loadfile,	   /* Optional, may be NULL 
			      Filename to load dotplot from */

	char *mtxfile,	   /* Optional, may be NULL 
			      Filename to load score matrix from */

	char *featurefile, /* Optional, may be NULL 
			      Filename to load features from */

	float memoryLimit, /* Optional, may be NULL 
			      Maximum Mb allowed for dotplot */

	int   zoomFac,	   /* Optional, may be NULL
			      Compression of dotplot {1, 2, 3, ... }
			      Automatically calculated if NULL */

	MSP  *MSPlist,	   /* Optional, may be NULL
			      List of MSPs containing genes and blast matches */

	int   MSPoffset,   /* Optional, may be NULL
			      Coordinate offset of MSPs */

	char *winsize,	   /* String determining the window size */

	int   pixelFacset  /* Preset pixel factor */
);
