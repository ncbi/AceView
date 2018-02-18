/*  File: fingerp.h
 *  Author: Ulrich Sauvage (ulrich@kaa.cnrs-mop.fr) 
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1993
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 19 05:36 1994 (mieg)
 * Created: Wed Oct 13 12:56:32 1993 (ulrich)
 *-------------------------------------------------------------------
 */

/* $Id: fingerp.h,v 1.1.1.1 2002/07/19 20:23:19 sienkiew Exp $ */

#ifndef DEFINE_fingerp
#define DEFINE_fingerp

  /* donnees is a boolean matrix clones/lengths */
  /* donnees must be created as Array(unsigned char) before call */
void pmapCptGetData (Array marInDef, Array defInMar, KEYSET clones, KEYSET bands) ;
void fpClearDisplay (void) ;
void fpDisplay (KEY clone) ;

#endif
