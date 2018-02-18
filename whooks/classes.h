/*  File: classes.h
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  1 11:42 1996 (srk)
 * Created: Thu Sep  1 15:59:43 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: classes.h,v 1.5 2006/10/06 03:23:17 mieg Exp $ */

  /* wspec/sysclass.h
     
	This file is read only at compile time

     It holds the enumeration of the classes needed by the applications.
     The classes, declared here as external, are explicitelly
     allocated and initialised in w4/tags.c 

     IF YOU EDIT this file EDIT ALSO w4/tags.c.

     Classes present in the models but not used explicitelly by the code
     no longer need to be declared.
     However, for backwards compatibility, code 3-* will still
     read the file wspec/classes.wrm if this file exist.

     If a class is needed in a single module, it does not have to be declared
     here but can as well be declared and initialised in that module
     following the example of tags.c
  */


extern int 
  _VPeptide,
  _VSequence,
  _VProtein,
  _VDNA,
  _VPaper,
  _VMethod,
  _VMap,
  _VgMap,
  _VvMap,
  _VMultiMap,
  _VLocus,
  _VGene,
  _VAllele,
  _VInterval,
  _V2_point_data,
  _VMulti_pt_data,
  _VClone,
  _VClone_Grid,
  _VPool,
  _VContig,
  _VpMap,
  _VChrom_Band,
  _VMotif,
  _VBaseCall,
  _VBaseQuality,
  _VBasePosition,
  _VOligoRepeat,
  _VHomology_group,
  _VGenetic_code,
  _VMap_set,
  _VDoc,
  _VPerson,
  _VPfam
;
 
 

