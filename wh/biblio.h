/*  File: biblio.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * $Id: biblio.h,v 1.1.1.1 2002/07/19 20:23:18 sienkiew Exp $
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 25 21:29 1996 (rd)
 * Created: Thu Jan 25 21:29:10 1996 (rd)
 *-------------------------------------------------------------------
 */

#ifndef DEFINE_BIBLIO_h
#define DEFINE_BIBLIO_h

  /* graphic func: display biblio associated to key or keyset */
void biblioKey(KEY key) ;
void biblioKeySet (char *title, KEYSET s) ;
  /* general: available from command.c */
void biblioDump (KEYSET bibSet, BOOL abstracts, int width) ; /* dumps on freeOut */
KEYSET biblioFollow (KEYSET s) ; /* associated papers etc */
BOOL biblioKeyPossible (KEY k) ; /* true if biblio tags in key */
BOOL biblioPossible (KEYSET s) ; /* true if biblio tags in class(keySet(s,0)) */

#endif
